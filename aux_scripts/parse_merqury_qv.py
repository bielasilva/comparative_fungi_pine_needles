#!/usr/bin/env python3
"""
parse_merqury_qv.py — aggregate per-isolate Merqury .qv files into Table S2.

Expects one <ISOLATE>.qv file per isolate, all in one flat directory (this
is where runMe_merqury.sh's FINAL_DIR copies them). Merqury's .qv format is
tab-separated, one row per assembly passed to merqury.sh (plus a "Both" row
if two assemblies were compared, which does not apply here) with columns:

    assembly_name  unique_kmers  common_kmers  QV  error_rate

- assembly_name : name merqury derived from the fasta filename
- unique_kmers  : k-mers found ONLY in the assembly, not in the read db
                  (Merqury's proxy for consensus errors)
- common_kmers  : k-mers found in both assembly and read db
- QV            : Phred-scaled consensus quality (higher = more accurate;
                  can be "inf" if unique_kmers is 0, i.e. no detected errors)
- error_rate    : per-base error rate implied by QV

Format confirmed against marbl/merqury source and the Galaxy Tool Shed
merqury wrapper's own test fixtures (iuc/merqury), not guessed.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

def parse_qv_file(path: Path) -> dict:
    """Read one isolate's .qv file. Returns the single assembly-level row.

    If merqury.sh was ever run with two assemblies (asm1 + asm2), it also
    emits a "Both" row; that shouldn't happen in this pipeline (one ONT-only
    assembly per isolate) but if it does, only the isolate-named row is
    used and a note flags the extra row so it isn't silently dropped.
    """
    lines = [l.strip() for l in path.read_text().splitlines() if l.strip()]
    if not lines:
        return {"error": "empty file"}

    rows = [l.split("\t") for l in lines]
    # Merqury's .qv has no header row - every line is data.
    data_rows = [r for r in rows if len(r) >= 5]

    if not data_rows:
        return {"error": f"no parseable rows (got {len(rows)} raw lines)"}

    isolate_id = path.stem
    matched = [r for r in data_rows if r[0] == isolate_id]
    chosen = matched[0] if matched else data_rows[0]

    extra_rows = len(data_rows) - 1
    note = "OK" if extra_rows == 0 else f"{extra_rows} extra row(s) in file, only first used"

    qv_raw = chosen[3]
    qv = float("inf") if qv_raw.lower() in ("inf", "+inf") else float(qv_raw)

    return {
        "unique_kmers": int(chosen[1]),
        "common_kmers": int(chosen[2]),
        "qv": qv,
        "error_rate": float(chosen[4]),
        "note": note,
    }


def main():

    qv_dir = Path("/home/aubgxs001/comparative_fungi/results/merqury")
    out_file = Path(f"{qv_dir}/Table_S2_QV.tsv")

    if not qv_dir.is_dir():
        sys.exit(f"ERROR: qv dir not found: {qv_dir}")

    qv_files = sorted(qv_dir.glob("*.qv"))
    if not qv_files:
        sys.exit(f"ERROR: no *.qv files found under {qv_dir}")

    rows = []
    for qv_file in qv_files:
        isolate_id = qv_file.stem
        parsed = parse_qv_file(qv_file)

        if "error" in parsed:
            rows.append({
                "isolate_id": isolate_id, "qv": None, "error_rate_pct": None,
                "unique_kmers": None, "common_kmers": None,
                "note": f"PARSE FAILED: {parsed['error']}",
            })
            continue

        rows.append({
            "isolate_id": isolate_id,
            "qv": parsed["qv"],
            "error_rate_pct": parsed["error_rate"] * 100,
            "unique_kmers": parsed["unique_kmers"],
            "common_kmers": parsed["common_kmers"],
            "note": parsed["note"],
        })

    fieldnames = ["isolate_id", "qv", "error_rate_pct", "unique_kmers",
                  "common_kmers", "note"]

    with open(out_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    finite_qv = [r["qv"] for r in rows if isinstance(r["qv"], float) and math.isfinite(r["qv"])]
    n_inf = sum(1 for r in rows if r["qv"] == float("inf"))
    n_failed = sum(1 for r in rows if r["note"] and r["note"].startswith("PARSE FAILED"))

    print(f"Wrote {len(rows)} isolates to {out_file}")
    if finite_qv:
        print(f"QV range (finite values only): {min(finite_qv):.1f}-{max(finite_qv):.1f} "
              f"(mean {sum(finite_qv)/len(finite_qv):.1f})  <- drop into the '[range]' "
              f"placeholder in the Methods paragraph")
    if n_inf:
        print(f"NOTE: {n_inf} isolate(s) returned QV=inf (zero unique/error k-mers detected). "
              f"This is a real Merqury output, not a bug, but it usually means either a very "
              f"clean assembly or a read set too shallow/small to catch errors - sanity-check "
              f"coverage for these before reporting inf in the table.", file=sys.stderr)
    if n_failed:
        print(f"WARNING: {n_failed} isolate(s) failed to parse - check their .qv file directly.",
              file=sys.stderr)


if __name__ == "__main__":
    main()