#!/usr/bin/env python3
"""
parse_antismash_json.py — per-genome BGC class counts from antiSMASH JSON
(chapter plan Fig 5: "Secondary metabolite cluster classes per genome").

WHY THIS PARSES THE 'region' FEATURE
-------------------------------------
Verified against a real antiSMASH JSON output (fungi_BS00_CNT_190-3.json,
21 detected regions), not assumed from documentation. Each detected BGC is
represented by three feature types on the record: 'protocluster',
'cand_cluster', and 'region', all sharing the same genomic coordinates and
in the sample inspected, a clean 1:1:1 correspondence (21 of each). 'region'
is antiSMASH's own FINAL classification after combining/deduplicating the
lower-level detections, and carries the 'product' qualifier directly, so it
is parsed rather than reconstructing the same information from the lower
two feature types.

'areas' and modules.antismash.detection.hmm_detection.rule_results were
both empty on every record inspected (including ones with real regions), so
they are not used as a data source here despite appearing in the schema.

HYBRID CLUSTERS
---------------
A region's 'product' qualifier is a LIST, and can legitimately contain more
than one entry (verified in the real data: ['NRPS', 'T1PKS'] and
['NRP-metallophore', 'NRPS']). This script counts BOTH ways and writes both
columns, rather than picking one silently:
    n_regions          - one row per region regardless of product count
                          (the region itself, e.g. "1 hybrid NRPS/T1PKS
                          cluster" counts as 1 toward NRPS and 1 toward
                          T1PKS here too, OR toward a combined "NRPS-T1PKS"
                          label if --hybrid-mode=joint)
    n_regions_per_1000  - normalized version of the above
Default is --hybrid-mode=split (count once per product, standard convention
in BGC survey papers, consistent with how CAZy multi-group families are
already handled elsewhere in this pipeline). --hybrid-mode=joint keeps
hybrids as their own combined category (e.g. "NRPS-T1PKS") instead of
splitting them, if you decide that is more informative for Fig 5.

CONTIG-EDGE FLAG
-----------------
Each region's 'contig_edge' qualifier ("True"/"False") is preserved in the
long-format output. A cluster flagged True is truncated by a contig
boundary - a real caveat for fragmented ONT-only assemblies, worth a
supplementary note or a sensitivity check (Fig 5 with/without edge-truncated
clusters) rather than silently treating a partial cluster the same as a
complete one.

Usage:
    # single genome
    python3 parse_antismash_json.py --json fungi_SAMPLE.json --sample SAMPLE \\
        --out-dir out/

    # all genomes at once (glob)
    python3 parse_antismash_json.py \\
        --json-glob '/path/to/antismash_json/fungi_*.json' \\
        --out-dir out/ --retained-samples files/retained_samples.txt
"""

import argparse
import json
import re
import sys
from collections import Counter
from glob import glob
from pathlib import Path

import pandas as pd

# antiSMASH product labels that should be treated as facets of the SAME BGC
# class rather than distinct classes for Fig 5. Left empty by default -
# fill in only if you decide e.g. 'NRPS' and 'NRPS-like' should be merged.
# Deliberately not pre-populated: collapsing categories is a Results-section
# decision, not something to bake into the parser silently.
PRODUCT_ALIASES = {
    # "NRPS-like": "NRPS",
}


def infer_sample_from_filename(path: Path) -> str:
    """fungi_<SAMPLE>.json -> <SAMPLE>. Falls back to the bare stem."""
    m = re.match(r"^fungi_(.+)$", path.stem)
    return m.group(1) if m else path.stem


def normalize_product(name: str) -> str:
    return PRODUCT_ALIASES.get(name, name)


def parse_one_json(path: Path, sample: str):
    """Returns a list of region dicts: sample, record_id, region_number,
    products (list[str]), contig_edge (bool), location (str).
    """
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError as e:
        print(f"ERROR: {path} is not valid JSON: {e}", file=sys.stderr)
        return []

    records = data.get("records", [])
    if not records:
        print(f"WARNING: {path} has no 'records' key or an empty list.", file=sys.stderr)
        return []

    rows = []
    for record in records:
        for feat in record.get("features", []):
            if feat.get("type") != "region":
                continue
            q = feat.get("qualifiers", {})
            products = [normalize_product(p) for p in q.get("product", [])]
            if not products:
                print(f"WARNING: {sample}/{record.get('id','?')} has a region feature "
                      f"with no 'product' qualifier: {q}", file=sys.stderr)
                continue
            contig_edge = q.get("contig_edge", ["False"])[0].strip().lower() == "true"
            region_number = q.get("region_number", ["?"])[0]
            rows.append({
                "sample": sample,
                "record_id": record.get("id", "?"),
                "region_number": region_number,
                "location": feat.get("location", ""),
                "products": products,
                "n_products": len(products),
                "is_hybrid": len(products) > 1,
                "contig_edge": contig_edge,
            })
    return rows


def load_retained_samples(path):
    if path is None:
        return None
    ids = {l.split("#", 1)[0].strip() for l in Path(path).read_text().splitlines()}
    ids.discard("")
    if not ids:
        raise SystemExit(f"ERROR: {path} contained no sample IDs.")
    return ids


def build_tables(region_rows, total_genes_by_sample, hybrid_mode):
    """From the flat region list, build:
      - regions_long: one row per region (for a supplementary table / QC)
      - class_long:   one row per (sample, bgc_class) count, for Fig 5
    """
    regions_long = pd.DataFrame(region_rows)
    if regions_long.empty:
        return regions_long, pd.DataFrame(
            columns=["sample", "bgc_class", "n_regions", "total_genes", "n_regions_per_1000"])

    regions_long["joint_label"] = regions_long["products"].apply(lambda p: "-".join(p))

    class_counter = Counter()
    for row in region_rows:
        sample = row["sample"]
        if hybrid_mode == "joint":
            class_counter[(sample, "-".join(row["products"]))] += 1
        else:  # split: one count per product in the region
            for p in row["products"]:
                class_counter[(sample, p)] += 1

    class_rows = []
    all_samples = sorted({s for s, _ in class_counter} | set(total_genes_by_sample))
    all_classes = sorted({c for _, c in class_counter})
    for sample in all_samples:
        total_genes = total_genes_by_sample.get(sample, 0)
        for bgc_class in all_classes:
            n = class_counter.get((sample, bgc_class), 0)
            class_rows.append({
                "sample": sample,
                "bgc_class": bgc_class,
                "n_regions": n,
                "total_genes": total_genes,
                "n_regions_per_1000": (n / total_genes * 1000) if total_genes else 0,
            })
    class_long = pd.DataFrame(class_rows)
    return regions_long, class_long


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--json", type=Path, help="A single antiSMASH JSON file")
    src.add_argument("--json-glob", help="Glob matching multiple antiSMASH JSON files")
    ap.add_argument("--sample", default=None,
                     help="Sample name for --json. Required with --json unless "
                          "the filename follows fungi_<SAMPLE>.json.")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--retained-samples", type=Path, default=None,
                     help="One sample ID per line. Same allowlist used by "
                          "parse_funannotate_gff3_annotations.py.")
    ap.add_argument("--gene-counts", type=Path, default=None,
                     help="Optional TSV/CSV with 'sample' and 'total_genes' columns "
                          "(e.g. product_summary.csv from the funannotate parser) "
                          "to compute per-1000-genes. Without it, only raw counts "
                          "are written.")
    ap.add_argument("--hybrid-mode", choices=["split", "joint"], default="split",
                     help="split (default): a hybrid region counts once per product "
                          "it contains. joint: hybrids form their own combined "
                          "category (e.g. 'NRPS-T1PKS').")
    args = ap.parse_args()

    if args.json:
        files = [(args.json, args.sample or infer_sample_from_filename(args.json))]
    else:
        files = [(Path(p), infer_sample_from_filename(Path(p))) for p in glob(args.json_glob)]
        if not files:
            raise SystemExit(f"ERROR: no files matched {args.json_glob!r}")

    retained = load_retained_samples(args.retained_samples)
    if retained is not None:
        before = len(files)
        files = [(p, s) for p, s in files if s in retained]
        print(f"Sample filter: {len(files)}/{before} JSON file(s) retained.")

    total_genes_by_sample = {}
    if args.gene_counts:
        sep = "\t" if str(args.gene_counts).endswith((".tsv", ".tab")) else ","
        gc = pd.read_csv(args.gene_counts, sep=sep)
        gene_col = next((c for c in ("total_genes", "total_mrna_features") if c in gc.columns),
                        None)
        if gene_col is None:
            raise SystemExit(f"ERROR: {args.gene_counts} has neither 'total_genes' nor "
                             f"'total_mrna_features'. Columns present: {list(gc.columns)}")
        total_genes_by_sample = dict(zip(gc["sample"], gc[gene_col]))
    else:
        print("NOTE: no --gene-counts given; only raw region counts will be written "
              "(no per-1000-genes normalization).", file=sys.stderr)

    all_rows = []
    for path, sample in files:
        if not path.exists():
            print(f"WARNING: {path} does not exist, skipping.", file=sys.stderr)
            continue
        rows = parse_one_json(path, sample)
        print(f"  {sample}: {len(rows)} region(s) ({sum(r['is_hybrid'] for r in rows)} hybrid)")
        all_rows.extend(rows)

    if not all_rows:
        raise SystemExit("ERROR: no regions parsed from any file. Check the JSON structure "
                         "matches what this script expects (region feature with a "
                         "'product' qualifier).")

    missing_genes = sorted({r["sample"] for r in all_rows} - set(total_genes_by_sample))
    if missing_genes and args.gene_counts:
        print(f"WARNING: {len(missing_genes)} sample(s) parsed from antiSMASH have no "
              f"matching gene count and will show 0 for per-1000-genes: {missing_genes}",
              file=sys.stderr)

    regions_long, class_long = build_tables(all_rows, total_genes_by_sample, args.hybrid_mode)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    regions_long.to_csv(args.out_dir / "antismash_regions_long.csv", index=False)
    class_long.to_csv(args.out_dir / "antismash_bgc_class_counts.csv", index=False)

    n_hybrid = sum(1 for r in all_rows if r["is_hybrid"])
    n_edge = sum(1 for r in all_rows if r["contig_edge"])
    print(f"\nParsed {len(all_rows)} region(s) across {len(set(r['sample'] for r in all_rows))} "
          f"sample(s).")
    print(f"  {n_hybrid} hybrid region(s) (hybrid-mode={args.hybrid_mode})")
    print(f"  {n_edge} region(s) flagged contig_edge=True (truncated by a contig boundary)")
    print(f"Wrote:\n  {args.out_dir / 'antismash_regions_long.csv'}\n"
          f"  {args.out_dir / 'antismash_bgc_class_counts.csv'}")


if __name__ == "__main__":
    main()