#!/usr/bin/env python3
"""
filter_secretome_dtm2.py — build the EffectorP input from SignalP +
DeepTMHMM2 results (chapter plan 2.8 / secretome-effector step).

Written against the ACTUAL DeepTMHMM2 output, not its README:
  * `segments` entries are 3-element ARRAYS [name, start, end], not objects
    with name/start/end fields as the README states.
  * The last element of predictions.json is a metadata object with no "id";
    it is skipped.
  * Coordinates are 1-based and inclusive.

WHY THIS PARSES JSON AND NOT THE TOPOLOGY STRING
------------------------------------------------
DeepTMHMM2's topology_string uses MEMBRANE-SPECIFIC characters by default
(its --simplify-io flag switches to generic I/O). A real run uses well over
30 distinct characters, including 'M', 'B', 'S', 'I' and 'O' for things that
are NOT what those letters meant in DeepTMHMM1. Scanning that string for
"M"/"B" silently mis-calls segments. The named segments are unambiguous.

SEGMENT VOCABULARY (verified against a real 839-protein run)
------------------------------------------------------------
  membrane-spanning : TMhelix, Beta sheet
  non-spanning      : interfacial, reentrant
  sorting signals   : signal, transit_peptide
  soluble sides     : outside, inside, Cytoplasmic, Lumenal, Extracellular,
                      Periplasmic, Mitochondrial matrix, Intermembrane space,
                      Nucleus matrix, Virion surface, Intravirion,
                      Thylakoid space, Stroma
Names are matched exactly (case-insensitively). Anything unrecognized is
reported loudly and treated as non-disqualifying - never silently guessed.
Exact matching is deliberate: substring matching would classify
"Intermembrane space" as a transmembrane segment because it contains
"membrane".

WHAT COUNTS AS DISQUALIFYING
----------------------------
EffectorP's secretome definition is signal peptide present AND no
transmembrane domain. Only membrane-SPANNING segments (TMhelix, Beta sheet)
count. Re-entrant regions and interfacial helices are membrane-associated
but do not traverse the bilayer, so they are NOT disqualifying by default;
--exclude-reentrant / --exclude-interfacial flip that.

MATURE vs FULL-LENGTH INPUT
---------------------------
If you ran DeepTMHMM2 on SignalP's *processed* FASTA (signal peptides
already removed), pass --mature-input. Every residue is then already mature,
so any TM segment is genuine and no cleavage-site logic is needed. Detected
automatically-ish via a warning if few sequences start with M.

If you ran it on FULL-LENGTH sequences, leave --mature-input off: TM
segments are then only counted when they extend past the signal peptide,
because a signal peptide misread as a TM helix would otherwise discard a
legitimately secreted protein.

Usage:
    # 1. inventory segment names actually present
    python3 filter_secretome_dtm2.py --predictions predictions.json --list-segment-names

    # 2a. mature (SP-stripped) input
    python3 filter_secretome_dtm2.py --predictions predictions.json \\
        --fasta mature_secretome.fasta --mature-input \\
        --out-fasta secretome_noTM_for_effectorp.fasta \\
        --out-table secretome_TM_filter_summary.tsv

    # 2b. full-length input
    python3 filter_secretome_dtm2.py --predictions predictions.json \\
        --fasta signalp_positive_fulllength.fasta \\
        --signalp-summary output_signalp.gff3 \\
        --out-fasta secretome_noTM_for_effectorp.fasta \\
        --out-table secretome_TM_filter_summary.tsv
"""

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

SPANNING = {"tmhelix", "beta sheet"}
REENTRANT = {"reentrant"}
INTERFACIAL = {"interfacial"}
SIGNAL = {"signal"}
TRANSIT = {"transit_peptide"}
SIDES = {
    "outside", "inside", "cytoplasmic", "lumenal", "extracellular",
    "periplasmic", "mitochondrial matrix", "intermembrane space",
    "nucleus matrix", "virion surface", "intravirion", "thylakoid space",
    "stroma", "mitochondrial intermembrane space",
}


def classify_segment(name: str) -> str:
    n = name.strip().lower()
    if n in SPANNING:
        return "spanning_tm"
    if n in REENTRANT:
        return "reentrant"
    if n in INTERFACIAL:
        return "interfacial"
    if n in SIGNAL:
        return "signal_peptide"
    if n in TRANSIT:
        return "transit_peptide"
    if n in SIDES:
        return "side"
    return "UNKNOWN"


def norm_segment(seg):
    """DeepTMHMM2 writes [name, start, end]; tolerate dict form too."""
    if isinstance(seg, (list, tuple)) and len(seg) >= 3:
        return str(seg[0]), int(seg[1]), int(seg[2])
    if isinstance(seg, dict):
        return str(seg.get("name", "")), int(seg["start"]), int(seg["end"])
    raise ValueError(f"unrecognized segment structure: {seg!r}")


def load_predictions(path: Path) -> dict:
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        data = data.get("predictions", data.get("results", list(data.values())))
    records, skipped = {}, 0
    for entry in data:
        if not isinstance(entry, dict) or "id" not in entry:
            skipped += 1  # metadata footer
            continue
        pid = str(entry["id"]).split()[0]
        segs = []
        for seg in entry.get("segments", []) or []:
            try:
                segs.append(norm_segment(seg))
            except (ValueError, TypeError, KeyError, IndexError):
                print(f"WARNING: skipping malformed segment in {pid}: {seg!r}", file=sys.stderr)
        records[pid] = {
            "type": entry.get("type", ""),
            "segments": segs,
            "membrane_types": entry.get("membrane_types", []) or [],
        }
    if skipped:
        print(f"(skipped {skipped} non-prediction record(s), e.g. the metadata footer)",
              file=sys.stderr)
    return records


def parse_fasta(path: Path) -> dict:
    records, header, seq_lines = {}, None, []

    def flush():
        if header is not None:
            records[header.split()[0]] = {"header": header, "seq": "".join(seq_lines)}

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            flush()
            header, seq_lines = line[1:], []
        else:
            seq_lines.append(line.strip())
    flush()
    return records


def parse_signalp_cleavage(path: Path) -> dict:
    cleavage = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 5 and parts[2].strip().lower() in ("signal_peptide", "signal peptide"):
            try:
                cleavage[parts[0].split()[0]] = int(parts[4]) + 1
                continue
            except ValueError:
                pass
        m = re.search(r"CS pos:\s*(\d+)-(\d+)", line)
        if m:
            pid = parts[0].split()[0] if parts and parts[0].strip() else line.split()[0]
            cleavage[pid] = int(m.group(2))
    return cleavage


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", required=True, type=Path)
    ap.add_argument("--list-segment-names", action="store_true")
    ap.add_argument("--fasta", type=Path,
                     help="FASTA whose surviving records are written out (the same one "
                          "you gave DeepTMHMM2, unless you intend to hand EffectorP "
                          "different sequences)")
    ap.add_argument("--out-fasta", type=Path)
    ap.add_argument("--out-table", type=Path)
    ap.add_argument("--mature-input", action="store_true",
                     help="Sequences already have signal peptides removed; count every "
                          "TM segment regardless of position")
    ap.add_argument("--signalp-summary", type=Path, default=None)
    ap.add_argument("--default-mature-start", type=int, default=25)
    ap.add_argument("--min-tm-overlap", type=int, default=5)
    ap.add_argument("--exclude-reentrant", action="store_true")
    ap.add_argument("--exclude-interfacial", action="store_true")
    args = ap.parse_args()

    preds = load_predictions(args.predictions)
    if not preds:
        sys.exit(f"ERROR: no predictions parsed from {args.predictions}")

    name_counter = Counter()
    for rec in preds.values():
        for name, _, _ in rec["segments"]:
            name_counter[name.strip()] += 1

    if args.list_segment_names:
        print(f"{len(preds)} predictions in {args.predictions}\n")
        print(f"{'count':>8}  {'classified as':<16}  segment name")
        for name, count in name_counter.most_common():
            print(f"{count:>8}  {classify_segment(name):<16}  {name!r}")
        print("\nprotein 'type' values:")
        for t, c in Counter(r["type"] for r in preds.values()).most_common():
            print(f"{c:>8}  {t!r}")
        unknown = [n for n in name_counter if classify_segment(n) == "UNKNOWN"]
        if unknown:
            print(f"\nWARNING: unrecognized segment name(s): {unknown}\nAdd them to the "
                  f"sets at the top of this script before filtering.", file=sys.stderr)
        return

    for req in ("fasta", "out_fasta", "out_table"):
        if getattr(args, req) is None:
            sys.exit(f"ERROR: --{req.replace('_','-')} required unless --list-segment-names")

    unknown = sorted(n for n in name_counter if classify_segment(n) == "UNKNOWN")
    if unknown:
        print(f"WARNING: unrecognized segment name(s) treated as NON-disqualifying: "
              f"{unknown}", file=sys.stderr)

    fasta_records = parse_fasta(args.fasta)
    cleavage = parse_signalp_cleavage(args.signalp_summary) if args.signalp_summary else {}

    # sanity check on mature vs full-length
    seqs = [r["seq"] for r in fasta_records.values()]
    frac_met = sum(1 for s in seqs if s.startswith("M")) / len(seqs) if seqs else 0
    if frac_met < 0.5 and not args.mature_input:
        print(f"WARNING: only {100*frac_met:.1f}% of sequences start with M, which suggests "
              f"signal peptides are already removed. If so, re-run with --mature-input; "
              f"otherwise TM segments near the N-terminus will be wrongly rescued.",
              file=sys.stderr)
    if frac_met > 0.9 and args.mature_input:
        print(f"WARNING: {100*frac_met:.1f}% of sequences start with M, which suggests these "
              f"are FULL-LENGTH sequences. --mature-input will then count signal peptides "
              f"misread as TM helices as real TM domains.", file=sys.stderr)

    pred_ids, fasta_ids = set(preds), set(fasta_records)
    if fasta_ids - pred_ids:
        m = sorted(fasta_ids - pred_ids)
        print(f"NOTE: {len(m)} FASTA ID(s) have no DeepTMHMM2 result; excluded. "
              f"First few: {m[:5]}", file=sys.stderr)
    if pred_ids - fasta_ids:
        e = sorted(pred_ids - fasta_ids)
        print(f"WARNING: {len(e)} DeepTMHMM2 ID(s) absent from --fasta (ID mismatch?). "
              f"First few: {e[:5]}", file=sys.stderr)

    rows, kept_ids, n_src = [], [], Counter()

    for pid in sorted(pred_ids & fasta_ids):
        rec = preds[pid]
        segs = [(classify_segment(n), n, s, e) for n, s, e in rec["segments"]]

        if args.mature_input:
            mature_start, source = 1, "mature_input"
        elif pid in cleavage:
            mature_start, source = cleavage[pid], "signalp_cs"
        else:
            sp = [x for x in segs if x[0] == "signal_peptide"]
            if sp:
                mature_start, source = max(x[3] for x in sp) + 1, "dtm2_signal_segment"
            else:
                mature_start, source = args.default_mature_start, "default_window"
        n_src[source] += 1

        spanning = [x for x in segs if x[0] == "spanning_tm"]
        spanning_mature = [x for x in spanning
                           if (x[3] - max(x[2], mature_start) + 1) >= args.min_tm_overlap]
        n_re = sum(1 for x in segs if x[0] == "reentrant")
        n_if = sum(1 for x in segs if x[0] == "interfacial")

        keep = not spanning_mature
        if args.exclude_reentrant and n_re:
            keep = False
        if args.exclude_interfacial and n_if:
            keep = False
        if keep:
            kept_ids.append(pid)

        rows.append({
            "protein_id": pid,
            "dtm2_type": rec["type"],
            "membrane_types": ";".join(map(str, rec["membrane_types"])),
            "mature_start": mature_start,
            "mature_start_source": source,
            "n_spanning_tm_total": len(spanning),
            "n_spanning_tm_counted": len(spanning_mature),
            "n_reentrant": n_re,
            "n_interfacial": n_if,
            "tm_rescued_as_signal_peptide": bool(spanning) and not spanning_mature,
            "spanning_tm_coords": ";".join(f"{x[2]}-{x[3]}" for x in spanning),
            "kept_for_effectorp": keep,
        })

    with open(args.out_fasta, "w") as fh:
        for pid in kept_ids:
            r = fasta_records[pid]
            fh.write(f">{r['header']}\n{r['seq']}\n")

    cols = ["protein_id", "dtm2_type", "membrane_types", "mature_start",
            "mature_start_source", "n_spanning_tm_total", "n_spanning_tm_counted",
            "n_reentrant", "n_interfacial", "tm_rescued_as_signal_peptide",
            "spanning_tm_coords", "kept_for_effectorp"]
    with open(args.out_table, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_total, n_kept = len(rows), len(kept_ids)
    if not n_total:
        sys.exit("ERROR: no proteins matched between the two files - check IDs.")

    print(f"{n_total} proteins with both a FASTA record and a DeepTMHMM2 result.")
    print(f"  {n_kept} kept ({100*n_kept/n_total:.1f}%) -> {args.out_fasta}")
    print(f"  {n_total - n_kept} dropped")
    print(f"  re-entrant in {sum(1 for r in rows if r['n_reentrant'])}, interfacial in "
          f"{sum(1 for r in rows if r['n_interfacial'])} protein(s) "
          f"[{'excluded' if args.exclude_reentrant else 'kept'} / "
          f"{'excluded' if args.exclude_interfacial else 'kept'}]")
    if not args.mature_input:
        print(f"  {sum(1 for r in rows if r['tm_rescued_as_signal_peptide'])} rescued "
              f"(TM call confined to the signal peptide)")
    print(f"  mature-start source: {dict(n_src)}")

    if n_kept / n_total < 0.5:
        print(f"\nNOTE: keeping under half is unusual for a fungal secretome (~80-90% "
              f"typical). Inspect spanning_tm_coords on dropped rows.", file=sys.stderr)


if __name__ == "__main__":
    main()
