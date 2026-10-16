#!/usr/bin/env python3
"""
Parse funannotate GFF3 annotation files and summarise functional content,
with a focus on CAZyme-based and EC-based substrate/activity grouping for a
comparative-genomics study of pine-needle Dothideomycetes.

This script MERGES two earlier scripts:
  - parse_funannotate_gff3_annotations.py  (rich annotation extraction:
        product, GO, PFAM, InterPro, EggNOG, COG, CAZy, plus per-sample
        summary tables and the output scaffolding / glob-based main())
  - parse_cazymes.py                       (EC-number extraction, the curated
        EC->activity and CAZy->activity maps, and the hybrid EC|CAZy evidence)

Per the request, CAZyme-related logic (normalisation, the CAZY_ACTIVITY_MAP,
family->group mapping) is taken from parse_cazymes.py.
"""

from pathlib import Path
from urllib.parse import unquote
from collections import Counter
from glob import glob
import argparse
import re
import sys

import pandas as pd


# ===========================================================================
# Category names (defined ONCE so the EC map and CAZy map cannot diverge).
# FIX 1: shared constants guarantee the hybrid union actually merges.
# ===========================================================================

CAT_CELLULOSE  = "Cellulose-associated (PCWDE)"
CAT_HEMICELL   = "Hemicellulose / xylan-associated (PCWDE)"
CAT_MANNAN     = "Mannan-associated (PCWDE)"
CAT_PECTIN     = "Pectin-associated (PCWDE)"
CAT_CUTIN      = "Cutin / suberin / cuticle esterase (infection-relevant)"
CAT_LPMO       = "Oxidative polysaccharide cleavage / LPMO"
CAT_REDOX      = "Auxiliary redox / lignin-related oxidation (mostly NOT lignin in Dothideomycetes)"
CAT_STARCH     = "Starch / storage glucan-associated"
CAT_CHITIN     = "Chitin / fungal cell wall - chitin (FCWE)"
CAT_GLUCAN_FCW = "Fungal cell wall - beta/alpha-glucan (FCWE, remodelling)"
UNGROUPED      = "Other / ungrouped CAZy"


# ===========================================================================
# COG category code -> description
# ===========================================================================

COG_MAP = {
    "J": "Translation ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control cell division chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall membrane envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking secretion and vesicular transport",
    "O": "Posttranslational modification protein turnover chaperones",
    "X": "Mobilome prophages transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


# ===========================================================================
# EC number -> activity group (high-specificity evidence).
# FIX 1: category names are the shared constants above.
# FIX 4: no greedy prefix map. Specific ester ECs replace the "3.1.1." blanket.
# ===========================================================================

EC_ACTIVITY_MAP = {
    # ---- Cellulose ----
    "3.2.1.4":   {CAT_CELLULOSE},   # endo-1,4-beta-glucanase (cellulase)
    "3.2.1.91":  {CAT_CELLULOSE},   # cellulose 1,4-beta-cellobiosidase, non-reducing end (CBH II-type)
    "3.2.1.176": {CAT_CELLULOSE},   # cellulose 1,4-beta-cellobiosidase, reducing end (CBH I-type)
    "3.2.1.21":  {CAT_CELLULOSE},   # beta-glucosidase (cellobiose->glucose); pairs with GH1/GH3.
                                    # POLYSPECIFIC in the same way those families are: broadens
                                    # coverage rather than sharpening specificity.
    # NOTE (FIX): 3.2.1.176 is a cellobiohydrolase, NOT an LPMO. Cellulose LPMOs
    # are oxidoreductases (EC 1.14.99.54 / .56) and are captured via CAZy AA9/AA16.

    # ---- Hemicellulose / xylan ----
    "3.2.1.8":   {CAT_HEMICELL},    # endo-1,4-beta-xylanase
    "3.2.1.37":  {CAT_HEMICELL},    # xylan 1,4-beta-xylosidase
    "3.2.1.72":  {CAT_HEMICELL},    # xylan 1,3-beta-xylosidase
    "3.2.1.55":  {CAT_HEMICELL},    # alpha-L-arabinofuranosidase (arabinoxylan; also pectic arabinan)
    "3.1.1.72":  {CAT_HEMICELL},    # acetylxylan esterase (specific; replaces "3.1.1." prefix)
    "3.1.1.73":  {CAT_HEMICELL},    # feruloyl esterase (arabinoxylan-lignin ester cross-links)

    # ---- Mannan ----
    "3.2.1.78":  {CAT_MANNAN},      # mannan endo-1,4-beta-mannosidase (beta-mannanase)
    "3.2.1.25":  {CAT_MANNAN},      # beta-mannosidase

    # ---- Pectin ----
    "3.2.1.15":  {CAT_PECTIN},      # endo-polygalacturonase
    "3.2.1.67":  {CAT_PECTIN},      # galacturan 1,4-alpha-galacturonidase (exo-PG)
    "3.2.1.99":  {CAT_PECTIN},      # arabinan endo-1,5-alpha-L-arabinanase  [MOVED from hemicellulose]
    "3.1.1.11":  {CAT_PECTIN},      # pectinesterase / pectin methylesterase (specific ester EC)
    "4.2.2.2":   {CAT_PECTIN},      # pectate lyase
    "4.2.2.10":  {CAT_PECTIN},      # pectin lyase

    # ---- Cutin / cuticle (infection-relevant) ----
    "3.1.1.74":  {CAT_CUTIN},       # cutinase

    # ---- Auxiliary redox / lignin-related (mostly NOT true ligninolysis here) ----
    "1.10.3.2":  {CAT_REDOX},       # laccase (AA1 multicopper oxidase)
    "1.11.1.13": {CAT_REDOX},       # manganese peroxidase (AA2)
    "1.11.1.14": {CAT_REDOX},       # lignin peroxidase (AA2)
    "1.11.1.16": {CAT_REDOX},       # versatile peroxidase (AA2)

    # ---- Starch / storage glucan ----
    "3.2.1.1":   {CAT_STARCH},      # alpha-amylase
    "3.2.1.2":   {CAT_STARCH},      # beta-amylase
    "3.2.1.3":   {CAT_STARCH},      # glucan 1,4-alpha-glucosidase (glucoamylase)
    "3.2.1.20":  {CAT_STARCH},      # alpha-glucosidase

    # ---- Chitin / fungal cell wall - chitin (FCWE) ----
    "3.2.1.14":  {CAT_CHITIN},      # chitinase
    "3.2.1.52":  {CAT_CHITIN},      # beta-N-acetylhexosaminidase

    # ---- Fungal cell wall - beta/alpha-glucan (FCWE) ----
    "3.2.1.6":   {CAT_GLUCAN_FCW},  # endo-1,3(4)-beta-glucanase  [MOVED out of mannan/beta-glucan]
    "3.2.1.39":  {CAT_GLUCAN_FCW},  # glucan endo-1,3-beta-D-glucosidase  [MOVED]
}

# Optional specific ECs you may add if your annotation uses them (kept out of the
# active map to avoid over-calling): 3.1.1.6 acetylesterase (too generic, skip),
# 3.5.1.41 chitin deacetylase -> CAT_CHITIN, 3.2.1.58 exo-1,3-beta-glucosidase ->
# CAT_GLUCAN_FCW, 2.4.1.25 4-alpha-glucanotransferase -> CAT_STARCH.


# ===========================================================================
# CAZy family -> activity group (taken from parse_cazymes.py, unchanged set).
# FIX 1: keys are the shared constants. Family lists and the biological comments
# are exactly as in the reviewed map.
# ===========================================================================

CAZY_ACTIVITY_MAP = {

    # ---------------- PLANT CELL WALL: structural polysaccharides ----------------

    CAT_CELLULOSE: {
        "GH5",    # endoglucanase; clan GH-A, HIGHLY POLYSPECIFIC -> use subfamily
        "GH6",    # cellobiohydrolase II (non-reducing end)
        "GH7",    # cellobiohydrolase I / endoglucanase (reducing end)
        "GH12",   # endoglucanase (also some xyloglucanase)
        "GH45",   # endoglucanase
        "GH48",   # cellobiohydrolase / reducing-end exo-glucanase (rare in fungi)
        "AA9",    # cellulose-active LPMO (formerly GH61)        # MULTI (LPMO)
        "CBM1",   # fungal cellulose-binding module               # MULTI (modular)
        "GH1",    # beta-glucosidase (cellobiose->glucose); POLYSPECIFIC
        "GH3",    # beta-glucosidase / beta-xylosidase; POLYSPECIFIC  # MULTI
    },

    CAT_HEMICELL: {
        "GH10", "GH11",          # endo-1,4-beta-xylanases
        "GH30",                  # glucuronoxylanase (subfamily-dependent)
        "GH51", "GH54",          # alpha-L-arabinofuranosidases   # MULTI (pectin side chains)
        "GH62",                  # arabinoxylan arabinofuranohydrolase
        "GH67", "GH115",         # alpha-glucuronidases
        "GH95",                  # alpha-L-fucosidase (xyloglucan/AGP)
        "GH27",                  # alpha-galactosidase (galacto(gluco)mannan side chains)
        "GH35",                  # beta-galactosidase             # MULTI (pectin galactan)
        "GH43",                  # beta-xylosidase / arabinofuranosidase; POLYSPECIFIC
        "GH74",                  # xyloglucanase (moved here from cellulose)
        "CE1",                   # acetyl xylan / feruloyl esterase; POLYSPECIFIC
        "CE2", "CE3", "CE6",     # acetyl xylan esterases
        "CE4",                   # acetyl xylan esterase          # MULTI (chitin deacetylase)
        "CE5",                   # acetyl xylan esterase          # MULTI (cutinase)
        "CE15",                  # glucuronoyl esterase (lignin-carbohydrate ester bonds)
        "CBM1",                  # cellulose-binding              # MULTI (modular)
    },

    CAT_MANNAN: {
        "GH26",                  # beta-mannanase                 # MULTI (hemicellulose)
        "GH5",                   # GH5_7 = beta-mannanase/mannosidase  # MULTI
        "GH2",                   # beta-mannosidase / beta-gal; POLYSPECIFIC  # MULTI
        "GH76",                  # alpha-1,6-mannanase
        "GH92", "GH38", "GH47",  # alpha-mannosidases
        "GH125",                 # exo-alpha-1,6-mannosidase
    },

    CAT_PECTIN: {
        "GH28",                  # polygalacturonase / rhamnogalacturonase
        "GH78",                  # alpha-L-rhamnosidase (RG-I)
        "GH88", "GH105",         # unsaturated (rhamno)galacturonyl hydrolases
        "GH53",                  # endo-beta-1,4-galactanase
        "GH93",                  # exo-arabinanase
        "GH79",                  # beta-glucuronidase
        "PL1", "PL3",            # pectate / pectin lyases
        "PL4", "PL9", "PL11",    # rhamnogalacturonan lyases (PL9 also pectate)
        "CE8",                   # pectin methylesterase
        "CE12",                  # pectin / rhamnogalacturonan acetylesterase  # MULTI
        "GH51", "GH54",          # arabinofuranosidases           # MULTI
        "GH35",                  # beta-galactosidase             # MULTI
    },

    # ---------------- CUTICLE / SUBERIN (lipid-polyester barriers) ----------------

    CAT_CUTIN: {
        "CE5",                   # fungal cutinase (EC 3.1.1.74) lives here  # MULTI
        # CAVEAT: CE5 also = acetyl xylan esterase; family signal cannot
        # distinguish cutinase vs AXE -> use phylogeny / subfamily / Pfam PF01083.
    },

    # ---------------- OXIDATIVE / AUXILIARY ACTIVITIES (redox) ----------------

    CAT_LPMO: {
        "AA9",                   # cellulose (fungi)              # MULTI
        "AA11",                  # chitin                         # MULTI (FCW)
        "AA13",                  # starch                         # MULTI (starch)
        "AA14",                  # xylan-coated cellulose (moved out of "lignin")
        "AA16",                  # cellulose
    },

    CAT_REDOX: {
        "AA1",   # multicopper oxidases / laccases; fungi usually AA1_3 (laccase-like)
        "AA2",   # class II peroxidases (LiP/MnP/VP) - LARGELY ABSENT in Dothideomycetes
        "AA3",   # GMC oxidoreductases (CDH, glucose oxidase, aryl-alcohol oxidase)
        "AA4",   # vanillyl-alcohol oxidase
        "AA5",   # copper radical oxidases (glyoxal oxidase) - H2O2 supply
        "AA6",   # 1,4-benzoquinone reductase (detox)
        "AA7",   # glucooligosaccharide oxidase
        "AA8",   # iron reductase domain (often fused to AA3)
    },

    # ---------------- STORAGE CARBOHYDRATE ----------------

    CAT_STARCH: {
        "GH13",                  # alpha-amylase
        "GH15",                  # glucoamylase
        "GH31",                  # alpha-glucosidase; POLYSPECIFIC  # MULTI (hemicell.)
        "GH77",                  # amylomaltase / 4-alpha-glucanotransferase
        "AA13",                  # starch-active LPMO              # MULTI
        "CBM20", "CBM21", "CBM25", "CBM26",
        "CBM34", "CBM41", "CBM48", "CBM53",
    },

    # ---------------- FUNGAL CELL WALL (self / antagonism, NOT plant biomass) ----------------

    CAT_CHITIN: {
        "GH18",                  # chitinase (endo-beta-N-acetylglucosaminidase)
        "GH19",                  # chitinase (mostly plant/bacterial; few fungi)
        "GH20",                  # beta-N-acetylhexosaminidase (NAGase)
        "GH75",                  # chitosanase
        "CE4",                   # chitin deacetylase             # MULTI (xylan AXE)
        "AA11",                  # chitin-active LPMO             # MULTI
        "CBM18",                 # chitin-binding module
        "CBM50",                 # LysM (chitin-binding; effector role, e.g. Ecp6)
    },

    CAT_GLUCAN_FCW: {
        "GH16",                  # clan GH-B; endo-1,3(4)-glucanase / XET; POLYSPECIFIC
        "GH17",                  # beta-1,3-glucanase / glucanosyltransferase
        "GH55",                  # exo/endo-beta-1,3-glucanase
        "GH64",                  # beta-1,3-glucanase
        "GH81",                  # endo-beta-1,3-glucanase
        "GH71",                  # alpha-1,3-glucanase (mutanase)
        "GH72",                  # beta-1,3-glucanosyltransglycosylase (cell-wall biogenesis)
    },
}


# Reverse index family -> {groups}, built once (efficient grouping).
FAMILY_TO_GROUPS = {}
for _group, _families in CAZY_ACTIVITY_MAP.items():
    for _fam in _families:
        FAMILY_TO_GROUPS.setdefault(_fam, set()).add(_group)


# FIX 1 (guard): EC labels must be a subset of the CAZy labels so they merge.
for _ec, _grps in EC_ACTIVITY_MAP.items():
    for _g in _grps:
        assert _g in CAZY_ACTIVITY_MAP, (
            f"EC {_ec} maps to category {_g!r} which is not a CAZy category. "
            f"EC and CAZy category names must match for the hybrid union."
        )


# ===========================================================================
# Attribute / value parsing helpers
# ===========================================================================

def parse_attributes(attr_string):
    """Parse a GFF3 column-9 attribute string into a dict.

    GFF3 attributes are ';'-separated key=value pairs; values may be URL-encoded
    and any literal ';' inside a value is percent-encoded as %3B, so splitting on
    ';' BEFORE unquoting is correct. Keys and items are stripped for safety.
    """
    attrs = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key.strip()] = unquote(value)
    return attrs


def get_attr(attrs, *keys, default=""):
    """Case-insensitive, first-match attribute lookup.

    FIX 6: protects against capitalisation differences ('note' vs 'Note',
    'Dbxref' vs 'dbxref') across funannotate versions / other annotators.
    """
    for key in keys:                       # exact match first (fast path)
        if key in attrs:
            return attrs[key]
    lower = {k.lower(): v for k, v in attrs.items()}
    for key in keys:
        if key.lower() in lower:
            return lower[key.lower()]
    return default


def split_csv_field(value):
    """Split a comma-separated field, dropping empties and '.' placeholders.

    Used for funannotate 'note' and 'Dbxref', which are comma-separated AFTER the
    attribute-level ';' split done by parse_attributes.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    value = str(value).strip()
    if value == "" or value == ".":
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def split_semicolon(value):
    """Split a ';'-joined field on ';' ONLY.

    FIX 2: activity-group names can contain commas
    (e.g. 'Fungal cell wall - beta/alpha-glucan (FCWE, remodelling)'), so group
    columns MUST be split on ';' only, never on ','.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    value = str(value).strip()
    if value == "" or value == ".":
        return []
    return [x.strip() for x in value.split(";") if x.strip()]


def normalize_cazy_family(fam):
    """Normalise a CAZy family string to family level.

    Examples: 'GH30_1' -> 'GH30', 'GH5.hmm' -> 'GH5', 'CAZy:GH7' -> 'GH7'.
    Subfamily detail is intentionally dropped here for activity grouping, but the
    RAW family string (with subfamily) is preserved in the output tables so it can
    be recovered (e.g. to separate GH5_7 mannanase from other GH5).
    """
    fam = str(fam).strip()
    fam = fam.replace("CAZy:", "")
    fam = fam.replace(".hmm", "")
    m = re.match(r"^(GH|GT|PL|CE|AA|CBM)(\d+)", fam)
    if m:
        return f"{m.group(1)}{m.group(2)}"
    return fam


def cazy_class(family):
    """Return the CAZy class prefix (GH/GT/PL/CE/AA/CBM) or 'Other'."""
    m = re.match(r"^(GH|GT|PL|CE|AA|CBM)", str(family))
    return m.group(1) if m else "Other"


def parse_dbxref(dbxref):
    """Split a Dbxref value into PFAM, InterPro, and other accessions."""
    pfams, interpros, others = [], [], []
    for item in split_csv_field(dbxref):
        if item.startswith("PFAM:"):
            pfams.append(item.replace("PFAM:", "", 1))
        elif item.startswith("InterPro:"):
            interpros.append(item.replace("InterPro:", "", 1))
        else:
            others.append(item)
    return pfams, interpros, others


def parse_note(note):
    """Parse a funannotate 'note' value, e.g. 'EggNog:ENOG503P00J,COG:G,CAZy:GH30'.

    Tag matching is case-insensitive (FIX 6). CAZy families are returned RAW
    (subfamily preserved); normalisation happens at point of use.
    """
    eggnog, cogs, cazys, other_notes = [], [], [], []
    for item in split_csv_field(note):
        up = item.upper()
        if up.startswith("EGGNOG:"):
            eggnog.append(item[len("EggNog:"):].strip())
        elif up.startswith("COG:"):
            cog_value = item[len("COG:"):].strip()
            for letter in cog_value:           # tolerate multi-letter codes
                if letter in COG_MAP:
                    cogs.append(letter)
        elif up.startswith("CAZY:"):
            cazy_value = item[len("CAZy:"):].strip()
            if cazy_value:
                cazys.append(cazy_value)
        else:
            other_notes.append(item)
    return eggnog, cogs, cazys, other_notes


def extract_ec_numbers(attrs):
    """Collect EC numbers from explicit EC keys and from note/Dbxref tokens.

    Handles EC_number / ec_number / EC keys (case-insensitive) and 'EC:' or
    'EC_number:' tokens embedded in note/Dbxref. Order-preserving de-duplication.
    """
    ecs = []
    for key, value in attrs.items():
        if key.lower() in {"ec_number", "ec"}:
            for v in split_csv_field(value):
                ecs.append(v.replace("EC:", "").replace("EC_number:", "").strip())
    for key, value in attrs.items():
        if key.lower() in {"note", "dbxref"}:
            for item in split_csv_field(value):
                for tag in ("EC_number:", "EC:"):
                    if item.startswith(tag):
                        ecs.append(item[len(tag):].strip())
                        break
    seen, out = set(), []
    for ec in ecs:
        if ec and ec not in seen:
            seen.add(ec)
            out.append(ec)
    return out


def ec_to_groups(ec):
    """Map a single EC number to activity groups (exact match only).

    FIX 4: no prefix matching. Incomplete ECs like '1.1.1.-' simply do not map,
    which is safer than catching catalase / lipase / glucose oxidase by accident.
    """
    return set(EC_ACTIVITY_MAP.get(ec, set()))


def cazy_to_groups(families, include_ungrouped=False):
    """Map an iterable of (normalised) CAZy families to a set of activity groups.

    Taken from parse_cazymes.py (interface preserved), backed by the reverse index
    for speed, with parse_funannotate's optional 'Other / ungrouped CAZy' bucket.
    """
    groups = set()
    for fam in families:
        g = FAMILY_TO_GROUPS.get(fam)
        if g:
            groups |= g
        elif include_ungrouped:
            groups.add(UNGROUPED)
    return groups


def infer_sample_name(gff_path):
    """Derive a sample name for one GFF3.

    FIX 7: funannotate names the annotate_results GFF3 <isolate>_<isolate>.gff3
    (e.g. BS00_120_204-2_BS00_120_204-2.gff3), so stripping only the extension
    yields a DOUBLED sample name ("BS00_120_204-2_BS00_120_204-2"). That name
    then fails to match the retained-sample allowlist and the secretome files,
    silently excluding every genome.

    Two defences:
      * from_directory=True takes the isolate directory name
        (.../funannotate_links/<SAMPLE>/predict/annotate_results/x.gff3),
        which is authoritative and immune to funannotate's file naming.
      * The filename fallback collapses an exactly-doubled stem.
    """

    name = gff_path.name
    name = re.sub(r"\.gff3?$", "", name)       # matches .gff and .gff3
    name = re.sub(r"\.annotations$", "", name)
    name = re.sub(r"\.final$", "", name)

    # collapse an exactly-doubled stem: "X_X" -> "X"
    half, rem = divmod(len(name) - 1, 2)
    if rem == 0 and name[half] == "_" and name[:half] == name[half + 1:]:
        name = name[:half]
    return name


# ===========================================================================
# Core GFF3 parsing
# ===========================================================================

HYPOTHETICAL_PRODUCTS = {"hypothetical protein", "hypothetical", "predicted protein"}


def parse_gff3(gff_path):
    """Parse one GFF3 file. Returns (list_of_mRNA_rows, feature_type_counts)."""
    rows = []
    feature_counts = Counter()
    sample = infer_sample_name(gff_path)

    with open(gff_path, "r", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attr_string = parts
            feature_counts[feature_type] += 1

            # mRNA carries the per-transcript functional annotation in funannotate.
            if feature_type != "mRNA":
                continue

            attrs = parse_attributes(attr_string)

            transcript_id = get_attr(attrs, "ID")
            gene_id = get_attr(attrs, "Parent")
            product = get_attr(attrs, "product", "Product")
            ontology_terms = split_csv_field(get_attr(attrs, "Ontology_term"))
            dbxref = get_attr(attrs, "Dbxref", "dbxref")
            note = get_attr(attrs, "note", "Note")

            pfams, interpros, other_dbxref = parse_dbxref(dbxref)
            eggnog, cogs, cazys_raw, other_notes = parse_note(note)
            cazys_norm = [normalize_cazy_family(c) for c in cazys_raw]

            ec_numbers = extract_ec_numbers(attrs)

            ec_groups = set()
            for ec in ec_numbers:
                ec_groups |= ec_to_groups(ec)
            cazy_groups = cazy_to_groups(cazys_norm)          # mapped groups only
            all_groups = ec_groups | cazy_groups              # FIX 1: now truly merges

            rows.append({
                "sample": sample,
                "gff_file": str(gff_path),
                "seqid": seqid,
                "source": source,
                "feature_type": feature_type,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "product": product,
                "is_hypothetical": product.lower() in HYPOTHETICAL_PRODUCTS,
                "go_terms": ";".join(ontology_terms),
                "pfams": ";".join(pfams),
                "interpros": ";".join(interpros),
                "eggnog_ids": ";".join(eggnog),
                "cog_categories": ";".join(cogs),
                "cazy_families": ";".join(cazys_raw),          # RAW (subfamily kept)
                "cazy_families_norm": ";".join(cazys_norm),    # family-level
                "ec_numbers": ";".join(ec_numbers),
                "ec_activity_groups": ";".join(sorted(ec_groups)),
                "cazy_activity_groups": ";".join(sorted(cazy_groups)),
                "all_activity_groups": ";".join(sorted(all_groups)),
                "has_ec_number": len(ec_numbers) > 0,
                "has_cazy": len(cazys_raw) > 0,
                "other_dbxref": ";".join(other_dbxref),
                "other_notes": ";".join(other_notes),
            })

    return rows, feature_counts


# ===========================================================================
# Summaries
# ===========================================================================

def count_multivalue(flat_df, column, sample_col="sample"):
    """Count ';'-separated occurrences of a column's values, per sample.

    This is OCCURRENCE counting (a gene with two PFAMs contributes to both). Used
    for PFAM / InterPro / GO / EggNOG / CAZy-family / COG, where that is intended.
    """
    rows = []
    for sample, sub in flat_df.groupby(sample_col):
        total_genes = sub.shape[0]
        counts = Counter()
        for value in sub[column].fillna(""):
            for item in split_semicolon(value):
                counts[item] += 1
        for item, count in counts.items():
            rows.append({
                "sample": sample,
                "feature": item,
                "count": count,
                "total_genes": total_genes,
                "count_per_1000_genes": count / total_genes * 1000 if total_genes else 0,
            })
    return pd.DataFrame(rows)


def activity_long(flat_df, group_col, evidence_label):
    """Explode a gene-level activity-group column into one row per (gene, group)."""
    rows = []
    for r in flat_df.itertuples(index=False):
        for group in split_semicolon(getattr(r, group_col)):
            rows.append({
                "sample": r.sample,
                "gene_id": r.gene_id,
                "transcript_id": r.transcript_id,
                "product": r.product,
                "activity_group": group,
                "evidence": evidence_label,
                "ec_numbers": r.ec_numbers,
                "cazy_families": r.cazy_families,
            })
    return pd.DataFrame(rows)


def activity_counts(long_df, total_genes_by_sample):
    """Count GENES per activity group, per sample.

    FIX 5: de-duplicate to (sample, transcript_id, activity_group) so a protein
    with two families in the same group is counted once. Column is 'n_genes'.
    """
    if long_df.empty:
        return pd.DataFrame(columns=[
            "sample", "activity_group", "n_genes", "total_genes", "n_genes_per_1000"
        ])
    counts = (
        long_df
        .drop_duplicates(["sample", "transcript_id", "activity_group"])
        .groupby(["sample", "activity_group"], as_index=False)
        .size()
        .rename(columns={"size": "n_genes"})
    )
    counts["total_genes"] = counts["sample"].map(total_genes_by_sample)
    counts["n_genes_per_1000"] = counts["n_genes"] / counts["total_genes"] * 1000
    return counts.sort_values(["sample", "activity_group"]).reset_index(drop=True)


def cazy_family_to_activity_long(cazy_family_counts):
    """Transparency table: per (sample, family) -> activity group(s), with the
    family OCCURRENCE count. Shows which families feed which groups. Not a gene
    count (use activity_counts for gene counts)."""
    rows = []
    for r in cazy_family_counts.itertuples(index=False):
        fam_norm = normalize_cazy_family(r.cazy_family)
        groups = FAMILY_TO_GROUPS.get(fam_norm, {UNGROUPED})
        for group in sorted(groups):
            rows.append({
                "sample": r.sample,
                "cazy_family": r.cazy_family,
                "cazy_family_norm": fam_norm,
                "activity_group": group,
                "family_occurrence_count": r.count,
                "total_genes": r.total_genes,
                "occ_per_1000_genes": r.count / r.total_genes * 1000 if r.total_genes else 0,
            })
    return pd.DataFrame(rows)


# ===========================================================================
# Sample filtering  (Methods 2.6: >15% BUSCO duplication exclusion)
# ===========================================================================

def load_retained_samples(path):
    """Read a one-ID-per-line allowlist. Blank lines and '#' comments ignored.

    Returns None when path is None, meaning 'no filtering'.
    """
    if path is None:
        return None
    ids = set()
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            ids.add(line)
    if not ids:
        raise SystemExit(f"ERROR: {path} contained no sample IDs.")
    return ids


def report_sample_filter(all_samples, retained):
    """Print kept/dropped samples and return the kept set.

    Prints rather than filtering silently, because the excluded set has to be
    stated in Methods and a typo in the allowlist would otherwise drop genomes
    invisibly.
    """
    if retained is None:
        print(f"No --retained-samples given: processing all {len(all_samples)} sample(s).")
        return set(all_samples)

    kept = {s for s in all_samples if s in retained}
    dropped = sorted(set(all_samples) - kept)
    missing = sorted(retained - set(all_samples))

    print(f"Sample filter: {len(kept)} kept, {len(dropped)} excluded.")
    if dropped:
        print(f"  excluded (not in allowlist): {dropped}")
    if missing:
        print(f"  WARNING: {len(missing)} allowlist ID(s) matched no GFF3: {missing}",
              file=sys.stderr)
    if not kept:
        raise SystemExit(
            "ERROR: sample filter removed every sample. Check that the allowlist ID "
            "format matches infer_sample_name() output (e.g. 'BS00_120_204-2')."
        )
    return kept


# ===========================================================================
# Activity grouping: inclusive (existing) vs exclusive (unambiguous only)
# ===========================================================================

FAMILY_GROUP_MULTIPLICITY = {fam: len(groups) for fam, groups in FAMILY_TO_GROUPS.items()}
MULTIGROUP_FAMILIES = {fam for fam, n in FAMILY_GROUP_MULTIPLICITY.items() if n > 1}


def annotate_family_table(cazy_family_df):
    """Add n_activity_groups / is_multigroup / is_ungrouped to the family table."""
    df = cazy_family_df.copy()
    norm = df["cazy_family"].apply(normalize_cazy_family)
    df["cazy_family_norm"] = norm
    df["n_activity_groups"] = norm.apply(lambda f: len(FAMILY_TO_GROUPS.get(f, ())))
    df["is_multigroup"] = df["n_activity_groups"] > 1
    df["is_ungrouped"] = df["n_activity_groups"] == 0
    return df


def exclusive_activity_long(flat_df):
    """Gene-level activity assignment counting ONLY unambiguous CAZy families.

    A family mapping to more than one activity group (CBM1 -> cellulose +
    hemicellulose, AA9 -> cellulose + LPMO, ...) is skipped rather than
    assigned to a guessed primary group. A gene is still counted if any of its
    OTHER families is unambiguous, so this removes ambiguous evidence, not
    whole genes.

    Anything present in both this table and the inclusive one is therefore not
    an artefact of multi-group membership.

    NOTE: reads the ';'-joined 'cazy_families_norm' column with
    split_semicolon. Using split_csv_field here would return the whole joined
    string as a single bogus family, and normalize_cazy_family would silently
    reduce "GH7;CBM1" to "GH7" -- a quiet, plausible-looking wrong answer.
    """
    rows = []
    for r in flat_df.itertuples(index=False):
        groups = set()
        for fam in split_semicolon(getattr(r, "cazy_families_norm", "")):
            fam_groups = FAMILY_TO_GROUPS.get(normalize_cazy_family(fam), set())
            if len(fam_groups) == 1:
                groups |= fam_groups
        for group in sorted(groups):
            rows.append({
                "sample": r.sample,
                "gene_id": r.gene_id,
                "transcript_id": r.transcript_id,
                "product": r.product,
                "activity_group": group,
                "evidence": "CAZy_family_unambiguous",
                "ec_numbers": r.ec_numbers,
                "cazy_families": r.cazy_families,
            })
    return pd.DataFrame(rows)


# ===========================================================================
# Secretome join  (filter_secretome_dtm2.py output)
# ===========================================================================

def substitute_locus_prefix(protein_id, prefix="FUN"):
    """Convert a funannotate locus-tag protein ID to its GFF3 transcript ID.

    The two numbering schemes are identical; only the prefix differs:
        BS00120204-2_000269-T1  ->  FUN_000269-T1

    Splits on the LAST underscore, so sample prefixes containing underscores
    are tolerated as long as the trailing '<number>-T<n>' part has none.
    """
    pid = str(protein_id).strip()
    if "_" not in pid:
        return pid
    return f"{prefix}_{pid.rsplit('_', 1)[-1]}"


ID_STRATEGIES = {
    "exact": lambda pid, prefix: str(pid).strip(),
    "locus_prefix": substitute_locus_prefix,
}


def load_secretome_dir(secretome_dir, pattern="*_secretome_TM_filter_summary.tsv"):
    """Load per-sample filter_secretome_dtm2.py summary tables.

    Sample ID is the filename with the pattern's fixed suffix removed.
    Returns {sample: {protein_id: kept_for_effectorp_bool}}.
    """
    secretome_dir = Path(secretome_dir)
    files = sorted(secretome_dir.glob(pattern))
    if not files:
        raise SystemExit(f"ERROR: no files matching {pattern!r} in {secretome_dir}")

    suffix = pattern.lstrip("*")
    out = {}
    for f in files:
        sample = f.name[:-len(suffix)] if f.name.endswith(suffix) else f.stem
        df = pd.read_csv(f, sep="\t")
        for col in ("protein_id", "kept_for_effectorp"):
            if col not in df.columns:
                raise SystemExit(f"ERROR: {f} lacks required column {col!r}")
        kept = df["kept_for_effectorp"].astype(str).str.strip().str.lower() == "true"
        out[sample] = dict(zip(df["protein_id"].astype(str).str.strip(), kept))
    print(f"Loaded secretome results for {len(out)} sample(s) from {secretome_dir}")
    return out


def load_id_map(path):
    """Optional override: two-column TSV, secretome_protein_id <TAB> gff_transcript_id."""
    if path is None:
        return None
    mapping = {}
    for i, line in enumerate(Path(path).read_text().splitlines()):
        if not line.strip() or line.startswith("#"):
            continue
        parts = [x.strip() for x in line.split("\t")]
        if len(parts) < 2:
            continue
        if i == 0 and parts[0].lower() in ("protein_id", "secretome_id"):
            continue
        mapping[parts[0]] = parts[1]
    if not mapping:
        raise SystemExit(f"ERROR: no usable pairs in {path}")
    return mapping


def attach_secretome(flat_df, secretome_by_sample, id_map=None, locus_tag_prefix="FUN"):
    """Add secretome_assessed / passed_tm_filter columns to the flat table.

    ID matching: for each sample every strategy in ID_STRATEGIES is tried
    against the transcript IDs actually present, and whichever matches most
    wins. This self-validates instead of trusting a naming rule that could
    differ per genome. An explicit id_map takes precedence when supplied.

    CRITICAL SEMANTICS: transcripts with no secretome record get pandas.NA,
    NOT False. A missing record means 'not assessed' (the protein was never in
    the SignalP-positive input), which is different from 'assessed and
    rejected'. Writing False would make every gene without a signal peptide
    look like a failed secretion candidate, and would silently zero the
    secreted-CAZyme counts if the join ever broke.
    """
    df = flat_df.copy()
    df["gff_id_tmp_"] = df["transcript_id"].astype(str).str.strip()
    gff_ids_by_sample = {s: set(sub["gff_id_tmp_"]) for s, sub in df.groupby("sample")}

    # Secretome files may exist for samples excluded by --retained-samples.
    # Ignore those rather than failing: they are not part of this analysis.
    extra = sorted(set(secretome_by_sample) - set(gff_ids_by_sample))
    if extra:
        print(f"  ignoring secretome file(s) for {len(extra)} sample(s) not in the "
              f"parsed set: {extra}")
    relevant = {s_: t for s_, t in secretome_by_sample.items() if s_ in gff_ids_by_sample}
    missing_sec = sorted(set(gff_ids_by_sample) - set(secretome_by_sample))
    if missing_sec:
        print(f"  WARNING: no secretome file for {len(missing_sec)} parsed sample(s): "
              f"{missing_sec}. Their CAZymes will all be 'not assessed'.",
              file=sys.stderr)

    remapped, chosen = {}, {}
    for sample, table in relevant.items():
        gff_ids = gff_ids_by_sample.get(sample, set())
        best_name, best_table, best_hits = None, {}, -1
        for name, fn in ID_STRATEGIES.items():
            candidate = {fn(pid, locus_tag_prefix): val for pid, val in table.items()}
            hits = len(set(candidate) & gff_ids)
            if hits > best_hits:
                best_name, best_table, best_hits = name, candidate, hits
        if id_map:
            forced = {id_map.get(pid, pid): val for pid, val in table.items()}
            hits = len(set(forced) & gff_ids)
            if hits >= best_hits:
                best_name, best_table, best_hits = "id_map", forced, hits
        remapped[sample] = best_table
        chosen[sample] = (best_name, best_hits, len(table))

    passed, assessed = [], []
    for r in df.itertuples(index=False):
        table = remapped.get(r.sample, {})
        if r.gff_id_tmp_ in table:
            passed.append(bool(table[r.gff_id_tmp_]))
            assessed.append(True)
        else:
            passed.append(pd.NA)
            assessed.append(False)

    df["secretome_assessed"] = assessed
    df["passed_tm_filter"] = passed
    df = df.drop(columns=["gff_id_tmp_"])

    print("Secretome join:")
    bad = []
    for sample, (strategy, hits, n_records) in sorted(chosen.items()):
        rate = hits / n_records if n_records else 0
        flag = ""
        if hits == 0:
            flag = "  <-- ZERO MATCHES"
            bad.append(sample)
        elif rate < 0.95:
            # Every secretome record should find its transcript, since those
            # proteins came from this annotation. A shortfall is ID drift,
            # not biology.
            flag = "  <-- some records unmatched"
        print(f"  {sample}: matched {hits}/{n_records} secretome records "
              f"({100*rate:.1f}%) via '{strategy}'{flag}")

    if bad:
        raise SystemExit(
            f"ERROR: {len(bad)} sample(s) had zero ID matches: {bad}\n"
            f"Neither exact matching nor '{locus_tag_prefix}_' prefix substitution linked "
            f"the secretome IDs to the GFF3 transcript IDs. Check --locus-tag-prefix, or "
            f"supply --secretome-id-map."
        )
    return df


def secreted_cazyme_summary(flat_df, total_genes_by_sample):
    """Per-sample CAZyme gene counts split by secretion status.

    The three status columns always sum to cazyme_genes.

    'not_assessed' means no predicted signal peptide, so the protein never
    entered the DeepTMHMM2 filter. 'nonsecreted' means it HAD a signal peptide
    but failed the TM filter, which is a much narrower claim. Report
    secreted/cazyme_genes, not secreted/(secreted+nonsecreted): the latter
    describes only the SignalP-positive subset and would badly overstate the
    secreted fraction.
    """
    rows = []
    for sample, sub in flat_df.groupby("sample"):
        total_genes = total_genes_by_sample.get(sample, 0)
        n_cazy = n_sec = n_nonsec = n_unassessed = 0
        for r in sub.itertuples(index=False):
            if not split_semicolon(getattr(r, "cazy_families", "")):
                continue
            n_cazy += 1
            if not getattr(r, "secretome_assessed", False):
                n_unassessed += 1
            elif bool(getattr(r, "passed_tm_filter")):
                n_sec += 1
            else:
                n_nonsec += 1
        rows.append({
            "sample": sample,
            "cazyme_genes": n_cazy,
            "secreted_cazyme_genes": n_sec,
            "nonsecreted_cazyme_genes": n_nonsec,
            "not_assessed_cazyme_genes": n_unassessed,
            "total_genes": total_genes,
            "cazy_per_1000_genes": n_cazy / total_genes * 1000 if total_genes else 0,
            "secreted_cazy_per_1000_genes": n_sec / total_genes * 1000 if total_genes else 0,
            "percent_of_cazymes_secreted": n_sec / n_cazy * 100 if n_cazy else 0,
        })
    return pd.DataFrame(rows)


def secreted_activity_counts(long_df, flat_df, total_genes_by_sample):
    """activity_counts() restricted to transcripts that passed the TM filter."""
    if long_df.empty:
        return pd.DataFrame(columns=["sample", "activity_group", "n_genes",
                                     "total_genes", "n_genes_per_1000"])
    secreted = flat_df.loc[
        flat_df["passed_tm_filter"].fillna(False).astype(bool),
        ["sample", "transcript_id"]
    ]
    keys = set(zip(secreted["sample"], secreted["transcript_id"]))
    mask = [ (s, t) in keys for s, t in zip(long_df["sample"], long_df["transcript_id"]) ]
    return activity_counts(long_df.loc[mask], total_genes_by_sample)


# ===========================================================================
# Wide matrices for the Fig 6 ordination
# ===========================================================================

def wide_matrix(long_df, columns, values, index="sample",
                normalize_by=None, fill=0):
    """Pivot a long count table into a sample x feature matrix.

    normalize_by: {sample: total_genes} converts to per-1000-genes. Use it.
    On raw counts the first PCA axis tends to track genome size rather than
    functional composition.
    """
    if long_df.empty:
        return pd.DataFrame()
    wide = long_df.pivot_table(index=index, columns=columns, values=values,
                               aggfunc="sum", fill_value=fill)
    if normalize_by:
        wide = wide.div(pd.Series(normalize_by).reindex(wide.index), axis=0) * 1000
    return wide.reset_index()


# ===========================================================================
# Main
# ===========================================================================


DEFAULT_GFF_DIR = ("/home/aubgxs001/scratch/comparative_fungi/results/funannote_parser/gff")
DEFAULT_SECRET_DIR = ("/home/aubgxs001/scratch/comparative_fungi/results/funannote_parser/dtm2")
DEFAULT_OUT_DIR = ("/home/aubgxs001/scratch/comparative_fungi/results/funannote_parser/results")
DEFAULT_SAMPLE_LIST = ("/home/aubgxs001/scratch/comparative_fungi/results/funannote_parser/samples_list.txt")


def build_parser():
    ap = argparse.ArgumentParser(
        description="Parse funannotate GFF3 annotations into CAZyme/EC activity tables.")
    ap.add_argument("--gff-dir", default=DEFAULT_GFF_DIR,
                    help="Directory with annotate_results GFF3 files (default: %(default)s)")
    ap.add_argument("--out-dir", type=Path, default=Path(DEFAULT_OUT_DIR))
    ap.add_argument("--retained-samples", type=Path, default=None,
                    help="One sample ID per line. Without it, ALL genomes are included. (default: %(default)s)")
    ap.add_argument("--secretome-dir", type=Path, default=DEFAULT_SECRET_DIR,
                    help="Directory of per-sample filter_secretome_dtm2.py summary TSVs. (default: %(default)s)")
    ap.add_argument("--secretome-pattern", default="*_secretome_TM_filter_summary.tsv")
    ap.add_argument("--secretome-id-map", type=Path, default=None,
                    help="Optional override TSV: secretome_protein_id <TAB> transcript_id")
    ap.add_argument("--locus-tag-prefix", default="FUN",
                    help="funannotate locus tag prefix used in the GFF3 (default: %(default)s)")
    ap.add_argument("--sample-from-filename", action="store_true",
                    help="Derive sample names from the GFF3 filename instead of the "
                         "isolate directory. Not recommended: funannotate writes "
                         "<isolate>_<isolate>.gff3, which doubles the name.")
    return ap


def main():
    args = build_parser().parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    gff_files = sorted(Path(args.gff_dir).glob("*.gff3"))

    # print(len(gff_files))
    # sys.exit("Stop")
    if not gff_files:
        raise SystemExit(
            f"ERROR: no GFF3 files matched:\n  {args.gff_dir}\n"
            f"Note that CAZy/COG/GO attributes exist only in annotate_results; a GFF3 "
            f"from predict_results parses fine but yields zero CAZymes."
        )

    from_dir = not args.sample_from_filename
    sample_of = {g: infer_sample_name(g) for g in gff_files}
    retained = load_retained_samples(args.retained_samples)
    keep = report_sample_filter(list(sample_of.values()), retained)
    gff_files = [g for g in gff_files if sample_of[g] in keep]

    all_rows = []
    feature_count_rows = []

    for gff in gff_files:
        sample = sample_of[gff]
        print(f"Parsing {sample}: {gff}")
        rows, feature_counts = parse_gff3(gff)
        for row in rows:                    # parse_gff3 re-derives from filename
            row["sample"] = sample
        all_rows.extend(rows)
        for feature, count in feature_counts.items():
            feature_count_rows.append({
                "sample": sample, "feature_type": feature,
                "count": count, "gff_file": str(gff),
            })

    flat = pd.DataFrame(all_rows)
    if flat.empty:
        raise SystemExit(
            "No mRNA features parsed. Check whether the GFF3 uses 'mRNA' features."
        )

    # ---- Item 6: secretome join (before the flat table is written) ----
    if args.secretome_dir:
        secretome = load_secretome_dir(args.secretome_dir, args.secretome_pattern)
        flat = attach_secretome(flat, secretome,
                                id_map=load_id_map(args.secretome_id_map),
                                locus_tag_prefix=args.locus_tag_prefix)
    else:
        flat["secretome_assessed"] = False
        flat["passed_tm_filter"] = pd.NA

    flat.to_csv(out_dir / "gene_annotation_flat_table.csv", index=False)
    pd.DataFrame(feature_count_rows).to_csv(
        out_dir / "gff3_feature_type_counts.csv", index=False
    )

    total_genes_by_sample = flat.groupby("sample").size().to_dict()

    # Item 3: the denominator is mRNA features. It equals the gene count only
    # when every gene has one transcript. Report any mismatch rather than
    # silently normalising by the wrong number.
    fc = pd.DataFrame(feature_count_rows)
    if not fc.empty and {"gene", "mRNA"} <= set(fc["feature_type"]):
        pivot = fc.pivot_table(index="sample", columns="feature_type",
                               values="count", aggfunc="sum", fill_value=0)
        if "gene" in pivot and "mRNA" in pivot:
            mismatch = pivot[pivot["gene"] != pivot["mRNA"]]
            if len(mismatch):
                print(f"WARNING: {len(mismatch)} sample(s) have mRNA != gene counts, so "
                      f"per-1000-genes values are per-1000-TRANSCRIPTS for those samples:",
                      file=sys.stderr)
                for smp, row in mismatch.iterrows():
                    print(f"    {smp}: {int(row['gene'])} genes, {int(row['mRNA'])} mRNAs",
                          file=sys.stderr)
            else:
                print("Denominator check: mRNA count == gene count for all samples.")

    # ---- Product / annotation-coverage summary (per sample) ----
    product_rows = []
    for sample, sub in flat.groupby("sample"):
        total = sub.shape[0]
        hypo = int(sub["is_hypothetical"].sum())
        product_rows.append({
            "sample": sample,
            "total_mrna_features": total,
            "annotated_product_count": total - hypo,
            "hypothetical_product_count": hypo,
            "hypothetical_percent": hypo / total * 100 if total else 0,
            "with_go_terms": int((sub["go_terms"] != "").sum()),
            "with_pfam": int((sub["pfams"] != "").sum()),
            "with_interpro": int((sub["interpros"] != "").sum()),
            "with_eggnog": int((sub["eggnog_ids"] != "").sum()),
            "with_cog": int((sub["cog_categories"] != "").sum()),
            "with_cazy": int(sub["has_cazy"].sum()),
            "with_ec": int(sub["has_ec_number"].sum()),
            "percent_with_cazy": int(sub["has_cazy"].sum()) / total * 100 if total else 0,
            "percent_with_ec": int(sub["has_ec_number"].sum()) / total * 100 if total else 0,
        })
    pd.DataFrame(product_rows).to_csv(out_dir / "product_summary.csv", index=False)

    # ---- COG category counts ----
    cog_rows = []
    for sample, sub in flat.groupby("sample"):
        total = sub.shape[0]
        counts = Counter()
        for value in sub["cog_categories"].fillna(""):
            for cog in split_semicolon(value):
                counts[cog] += 1
        for cog, count in counts.items():
            cog_rows.append({
                "sample": sample,
                "category_code": cog,
                "category": COG_MAP.get(cog, "Unknown COG category"),
                "count": count,
                "total_genes": total,
                "count_per_1000_genes": count / total * 1000 if total else 0,
            })
    pd.DataFrame(cog_rows).sort_values(["sample", "category_code"]).to_csv(
        out_dir / "cog_category_counts.csv", index=False
    ) if cog_rows else pd.DataFrame().to_csv(out_dir / "cog_category_counts.csv", index=False)

    # ---- CAZy family + class OCCURRENCE counts (raw family, subfamily kept) ----
    cazy_family = count_multivalue(flat, "cazy_families")
    if not cazy_family.empty:
        cazy_family = cazy_family.rename(columns={"feature": "cazy_family"})
        cazy_family["cazy_class"] = cazy_family["cazy_family"].apply(cazy_class)
        # Item 4: flag families that feed more than one activity group.
        cazy_family = annotate_family_table(cazy_family)
        cazy_family.to_csv(out_dir / "cazy_family_counts.csv", index=False)

        cazy_class_counts = (
            cazy_family.groupby(["sample", "cazy_class", "total_genes"], as_index=False)["count"].sum()
        )
        cazy_class_counts["count_per_1000_genes"] = (
            cazy_class_counts["count"] / cazy_class_counts["total_genes"] * 1000
        )
        cazy_class_counts.to_csv(out_dir / "cazy_class_counts.csv", index=False)

        # transparency: which families feed which activity groups (occurrence)
        cazy_family_to_activity_long(cazy_family).to_csv(
            out_dir / "cazy_family_to_activity_long.csv", index=False
        )
    else:
        for name in ["cazy_family_counts.csv", "cazy_class_counts.csv",
                     "cazy_family_to_activity_long.csv"]:
            pd.DataFrame().to_csv(out_dir / name, index=False)

    # ---- Activity-group GENE counts: CAZy, EC, and hybrid (FIX 5: gene-level) ----
    cazy_long = activity_long(flat, "cazy_activity_groups", "CAZy_family")
    ec_long = activity_long(flat, "ec_activity_groups", "EC")
    hybrid_long = activity_long(flat, "all_activity_groups", "EC_or_CAZy_family")

    cazy_long.to_csv(out_dir / "cazy_activity_gene_long.csv", index=False)
    ec_long.to_csv(out_dir / "ec_activity_gene_long.csv", index=False)
    hybrid_long.to_csv(out_dir / "hybrid_ec_cazy_activity_gene_long.csv", index=False)

    activity_counts(cazy_long, total_genes_by_sample).to_csv(
        out_dir / "cazy_activity_group_counts.csv", index=False)
    activity_counts(ec_long, total_genes_by_sample).to_csv(
        out_dir / "ec_activity_group_counts.csv", index=False)
    activity_counts(hybrid_long, total_genes_by_sample).to_csv(
        out_dir / "hybrid_ec_cazy_activity_group_counts.csv", index=False)

    # ---- Item 4: conservative companion using unambiguous families only ----
    excl_long = exclusive_activity_long(flat)
    excl_long.to_csv(out_dir / "cazy_activity_gene_long_exclusive.csv", index=False)
    excl_counts = activity_counts(excl_long, total_genes_by_sample)
    excl_counts.to_csv(out_dir / "cazy_activity_group_counts_exclusive.csv", index=False)

    # ---- Item 6: secreted subsets of the activity tables ----
    if args.secretome_dir:
        secreted_cazyme_summary(flat, total_genes_by_sample).to_csv(
            out_dir / "secreted_cazyme_summary.csv", index=False)
        secreted_activity_counts(cazy_long, flat, total_genes_by_sample).to_csv(
            out_dir / "cazy_activity_group_counts_secreted.csv", index=False)

    # ---- Wide matrices for the Fig 6 ordination ----
    if not cazy_family.empty:
        wide_matrix(cazy_family, columns="cazy_family_norm", values="count",
                    normalize_by=total_genes_by_sample).to_csv(
            out_dir / "matrix_cazy_family_per1000.csv", index=False)
        unamb = cazy_family[~cazy_family["is_multigroup"] & ~cazy_family["is_ungrouped"]]
        wide_matrix(unamb, columns="cazy_family_norm", values="count",
                    normalize_by=total_genes_by_sample).to_csv(
            out_dir / "matrix_cazy_family_unambiguous_per1000.csv", index=False)
    if not excl_counts.empty:
        wide_matrix(excl_counts, columns="activity_group", values="n_genes",
                    normalize_by=total_genes_by_sample).to_csv(
            out_dir / "matrix_activity_group_exclusive_per1000.csv", index=False)

    # ---- PFAM / InterPro / GO / EggNOG occurrence counts ----
    for column, outfile in [
        ("pfams", "pfam_counts.csv"),
        ("interpros", "interpro_counts.csv"),
        ("go_terms", "go_term_counts.csv"),
        ("eggnog_ids", "eggnog_counts.csv"),
    ]:
        counts = count_multivalue(flat, column)
        (counts if not counts.empty else pd.DataFrame()).to_csv(out_dir / outfile, index=False)

    print()
    print(f"Parsed {len(gff_files)} GFF3 file(s)")
    print(f"Wrote outputs to {out_dir}")
    print()
    print("Key outputs:")
    for f in [
        "gene_annotation_flat_table.csv",
        "product_summary.csv",
        "cog_category_counts.csv",
        "cazy_family_counts.csv",
        "cazy_class_counts.csv",
        "cazy_activity_group_counts.csv",
        "ec_activity_group_counts.csv",
        "hybrid_ec_cazy_activity_group_counts.csv",
        "cazy_activity_group_counts_exclusive.csv",
        "secreted_cazyme_summary.csv",
        "cazy_activity_group_counts_secreted.csv",
        "matrix_cazy_family_per1000.csv",
        "matrix_cazy_family_unambiguous_per1000.csv",
        "matrix_activity_group_exclusive_per1000.csv",
    ]:
        marker = "" if (out_dir / f).exists() else "   (not written)"
        print(f"  {out_dir / f}{marker}")


if __name__ == "__main__":
    main()