#!/usr/bin/env python3
"""
parse_effectorp.py — per-genome effector counts from EffectorP 3.0 TSV output.

Verified against real output (BS00_120_204-2.tsv, 708 proteins), not assumed
from the documented header.

TWO THINGS THE HEADER DOES NOT TELL YOU
---------------------------------------
1. The "# Identifier" column contains TWO space-separated IDs, transcript then
   gene:
       BS00120204-2_000021-T1 BS00120204-2_000021
   Splitting on tabs alone leaves them concatenated and nothing will join to
   the GFF3. Both are extracted here.

2. There are FIVE prediction values, not three. Dual localisation is written
   in both orders:
       Non-effector                     456
       Apoplastic effector              149
       Cytoplasmic effector              55
       Apoplastic/cytoplasmic effector   29
       Cytoplasmic/apoplastic effector   19
   The last two are the same biological call. They are normalised to a single
   "Dual" category, because splitting 48 genes across two buckets that mean
   the same thing would understate dual-localised effectors by whichever
   order happened to be less common.

PROBABILITY COLUMNS
-------------------
Each of the three probability columns holds either "-" or "Y (0.931)". The
numeric value is extracted where present. Note that a dual-prediction row
carries TWO probabilities (one cytoplasmic, one apoplastic), so there is no
single "the" effector probability for such genes; both are kept as separate
columns rather than collapsed.

ID CONVENTION
-------------
Transcript IDs use the funannotate locus-tag prefix
(BS00120204-2_000021-T1), while the annotate_results GFF3 uses FUN_ IDs
(FUN_000021-T1). Same numbering, different prefix. --locus-tag-prefix
handles the substitution, identically to filter_secretome_dtm2.py.

Usage:
    python3 parse_effectorp.py \\
        --effectorp-dir /path/to/effectorp/ \\
        --gene-counts files/product_summary.csv \\
        --retained-samples files/retained_samples.txt \\
        --out-dir files/
"""

import argparse
import re
import sys
from glob import glob
from pathlib import Path

import pandas as pd

# Normalised prediction categories. Both dual orderings collapse to "Dual".
PREDICTION_MAP = {
    "non-effector": "Non-effector",
    "apoplastic effector": "Apoplastic",
    "cytoplasmic effector": "Cytoplasmic",
    "apoplastic/cytoplasmic effector": "Dual",
    "cytoplasmic/apoplastic effector": "Dual",
}

PROB_RE = re.compile(r"Y\s*\(([0-9.]+)\)")


def parse_prob(cell):
    """'Y (0.931)' -> 0.931 ; '-' -> None. Anything else is reported."""
    cell = str(cell).strip()
    if cell in ("-", "", "nan"):
        return None
    m = PROB_RE.match(cell)
    if m:
        return float(m.group(1))
    print(f"WARNING: unrecognised probability cell {cell!r}", file=sys.stderr)
    return None


def substitute_locus_prefix(protein_id, prefix="FUN"):
    """BS00120204-2_000021-T1 -> FUN_000021-T1 (same numbering, new prefix)."""
    pid = str(protein_id).strip()
    return f"{prefix}_{pid.rsplit('_', 1)[-1]}" if "_" in pid else pid


def parse_one(path: Path, sample: str, locus_tag_prefix: str):
    rows = []
    unknown_preds = set()

    with open(path) as fh:
        header = fh.readline()
        if not header.startswith("#"):
            print(f"WARNING: {path.name} does not start with a '#' header line; "
                  f"first line may have been treated as data.", file=sys.stderr)
            fh.seek(0)

        for line_no, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                print(f"WARNING: {path.name}:{line_no} has {len(parts)} fields, "
                      f"expected 5. Skipping.", file=sys.stderr)
                continue

            # Identifier column: "<transcript_id> <gene_id>"
            ident = parts[0].strip().split()
            transcript_id = ident[0]
            gene_id = ident[1] if len(ident) > 1 else None

            pred_raw = parts[4].strip()
            pred = PREDICTION_MAP.get(pred_raw.lower())
            if pred is None:
                unknown_preds.add(pred_raw)
                pred = pred_raw   # keep verbatim rather than silently dropping

            rows.append({
                "sample": sample,
                "transcript_id": transcript_id,
                "gff_transcript_id": substitute_locus_prefix(transcript_id, locus_tag_prefix),
                "gene_id": gene_id,
                "prob_cytoplasmic": parse_prob(parts[1]),
                "prob_apoplastic": parse_prob(parts[2]),
                "prob_non_effector": parse_prob(parts[3]),
                "prediction_raw": pred_raw,
                "prediction": pred,
                "is_effector": pred != "Non-effector",
            })

    if unknown_preds:
        print(f"WARNING: {path.name} contains prediction value(s) not in "
              f"PREDICTION_MAP: {sorted(unknown_preds)}. Kept verbatim; add them "
              f"to the map before trusting the counts.", file=sys.stderr)
    return rows


def load_retained(path):
    if path is None:
        return None
    ids = {l.split("#", 1)[0].strip() for l in Path(path).read_text().splitlines()}
    ids.discard("")
    return ids or None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--effectorp-dir", type=Path, required=True)
    ap.add_argument("--pattern", default="*.tsv")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--retained-samples", type=Path, default=None)
    ap.add_argument("--gene-counts", type=Path, default=None,
                    help="CSV/TSV with 'sample' and 'total_genes' or "
                         "'total_mrna_features'. Enables per-1000 normalization.")
    ap.add_argument("--locus-tag-prefix", default="FUN")
    args = ap.parse_args()

    files = sorted(Path(p) for p in glob(str(args.effectorp_dir / args.pattern)))
    if not files:
        raise SystemExit(f"ERROR: no files matching {args.pattern!r} in {args.effectorp_dir}")

    retained = load_retained(args.retained_samples)
    if retained is not None:
        before = len(files)
        files = [f for f in files if f.stem in retained]
        print(f"Sample filter: {len(files)}/{before} file(s) retained.")
        missing = sorted(retained - {f.stem for f in files})
        if missing:
            print(f"  NOTE: {len(missing)} retained sample(s) have no EffectorP file: "
                  f"{missing[:5]}{'...' if len(missing) > 5 else ''}", file=sys.stderr)

    all_rows = []
    for f in files:
        rows = parse_one(f, f.stem, args.locus_tag_prefix)
        n_eff = sum(r["is_effector"] for r in rows)
        print(f"  {f.stem}: {len(rows)} proteins, {n_eff} effectors")
        all_rows.extend(rows)

    if not all_rows:
        raise SystemExit("ERROR: no rows parsed.")

    long_df = pd.DataFrame(all_rows)

    # --- per-genome summary ------------------------------------------------
    total_genes = {}
    if args.gene_counts:
        sep = "\t" if str(args.gene_counts).endswith((".tsv", ".tab")) else ","
        gc = pd.read_csv(args.gene_counts, sep=sep)
        col = next((c for c in ("total_genes", "total_mrna_features") if c in gc.columns), None)
        if col is None:
            raise SystemExit(f"ERROR: {args.gene_counts} has neither 'total_genes' nor "
                             f"'total_mrna_features'. Columns: {list(gc.columns)}")
        total_genes = dict(zip(gc["sample"], gc[col]))
    else:
        print("NOTE: no --gene-counts; per-1000 columns will be omitted.", file=sys.stderr)

    cats = ["Non-effector", "Apoplastic", "Cytoplasmic", "Dual"]
    summary = []
    for sample, sub in long_df.groupby("sample"):
        counts = sub["prediction"].value_counts().to_dict()
        n_screened = len(sub)
        n_eff = int(sub["is_effector"].sum())
        tg = total_genes.get(sample)
        row = {
            "sample": sample,
            "n_proteins_screened": n_screened,
            "n_effectors": n_eff,
            "n_apoplastic": counts.get("Apoplastic", 0),
            "n_cytoplasmic": counts.get("Cytoplasmic", 0),
            "n_dual": counts.get("Dual", 0),
            "n_non_effector": counts.get("Non-effector", 0),
            "pct_of_screened_effector": 100 * n_eff / n_screened if n_screened else 0,
            "total_genes": tg,
        }
        if tg:
            row["effectors_per_1000_genes"] = n_eff / tg * 1000
            row["apoplastic_per_1000_genes"] = counts.get("Apoplastic", 0) / tg * 1000
            row["cytoplasmic_per_1000_genes"] = counts.get("Cytoplasmic", 0) / tg * 1000
        summary.append(row)
    summary_df = pd.DataFrame(summary).sort_values("sample")

    # --- long-format category counts, for figures --------------------------
    cat_rows = []
    for sample, sub in long_df.groupby("sample"):
        tg = total_genes.get(sample)
        counts = sub["prediction"].value_counts().to_dict()
        for cat in cats:
            n = counts.get(cat, 0)
            cat_rows.append({
                "sample": sample, "effector_class": cat, "n_genes": n,
                "total_genes": tg,
                "n_genes_per_1000": (n / tg * 1000) if tg else None,
            })
    cat_df = pd.DataFrame(cat_rows)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(args.out_dir / "effectorp_predictions_long.csv", index=False)
    summary_df.to_csv(args.out_dir / "effectorp_summary.csv", index=False)
    cat_df.to_csv(args.out_dir / "effectorp_class_counts.csv", index=False)

    print(f"\nParsed {len(long_df)} proteins across {long_df['sample'].nunique()} sample(s).")
    print(f"  effectors: {int(long_df['is_effector'].sum())} "
          f"({100*long_df['is_effector'].mean():.1f}% of screened)")
    print(long_df["prediction"].value_counts().to_string())
    print(f"\nWrote:\n  {args.out_dir / 'effectorp_predictions_long.csv'}"
          f"\n  {args.out_dir / 'effectorp_summary.csv'}"
          f"\n  {args.out_dir / 'effectorp_class_counts.csv'}")

    if not total_genes:
        print("\nNOTE: without --gene-counts the only available denominator is "
              "n_proteins_screened, which is the SignalP-positive secretome, not the "
              "proteome. pct_of_screened_effector is therefore a fraction of the "
              "secretome. Supply --gene-counts for per-1000-genes values comparable "
              "to the CAZyme and BGC tables.", file=sys.stderr)


if __name__ == "__main__":
    main()
