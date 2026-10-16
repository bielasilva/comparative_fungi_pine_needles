# Comparative genomics of pine-needle fungal endophytes

Code for the paper "Functional gene content of fungi in a pine needle
disease system largely reflects their phylogeny"

A comparative genomic analysis of fungal endophytes isolated
from pine needles, alongside reference genomes of pine and plant pathogens,
saprophytes and yeasts. The analysis asks whether the functional repertoire
of each genome (CAZymes, secreted proteins and effectors, COG categories and
biosynthetic gene clusters) is structured by lifestyle or by taxonomy.

The repository contains the scripts that take the isolates from raw Nanopore
reads to annotated genomes, the equivalent annotation steps for the reference
genomes, the phylogenomic analysis, and the R/Quarto document with all of the
statistics and figures.

## Repository layout

```
.
├── isolates/                      Genomes sequenced in this study
│   ├── 01_assembly/               Flye -> Medaka (two rounds) -> LINKS scaffolding
│   ├── 02_assembly_qc/            QUAST, BUSCO, Merqury
│   ├── 03_repeat_masking/         RepeatMasker (soft-masking)
│   ├── 04_gene_prediction/        funannotate clean / sort / predict
│   ├── 05_secretome/              SignalP 6, DeepTMHMM2, EffectorP 3.0
│   ├── 06_functional_annotation/  InterProScan, antiSMASH, funannotate annotate
│   ├── 07_genome_similarity/      skani
│   └── utils/                     Run-all wrapper and progress counter
├── references/                    Public reference genomes (same step numbers as isolates/)
│   ├── 02_assembly_qc/            QUAST, BUSCO
│   ├── 04_gene_prediction/        funannotate clean / predict
│   ├── 05_secretome/              SignalP 6, DeepTMHMM2, EffectorP 3.0
│   └── 06_functional_annotation/  InterProScan, antiSMASH, funannotate annotate, OrthoFinder
├── phylogenomics/                 BUSCO single-copy genes -> supermatrix -> IQ-TREE
├── aux_scripts/                   Python parsers that turn tool outputs into tables
└── R/
    └── comparative_fungi_code.qmd   QC, taxonomic placement, statistics and figures
```

Step numbers follow the order in which the steps were run. The secretome step
(05) comes before functional annotation (06) because `funannotate annotate`
uses the SignalP and antiSMASH outputs. In `references/` the step numbers
match those in `isolates/`, so the gaps (no 01 or 03) are intentional: the
reference assemblies were downloaded, not assembled.

## Pipeline overview

### Isolates

| Step | Tool | Script |
|---|---|---|
| Assembly | Flye (`--nano-raw`) | `01_assembly/runMe_flye.sh` |
| Polishing | Medaka, two rounds (model `r1041_e82_400bps_sup_v5.2.0`) | `01_assembly/runMe_medaka01.sh`, `runMe_medaka02.sh` |
| Scaffolding | LINKS | `01_assembly/runMe_links.sh` |
| Assembly QC | QUAST, BUSCO (`fungi_odb12`), Merqury | `02_assembly_qc/` |
| Repeat masking | RepeatMasker (`-species fungi -xsmall`) | `03_repeat_masking/runMe_repeatmasker.sh` |
| Gene prediction | funannotate predict (`--busco_db dikarya`) | `04_gene_prediction/runMe_funannotate_links.sh` |
| Signal peptides | SignalP 6 | `05_secretome/runMe_signalp.sh` |
| Transmembrane filter | DeepTMHMM2 + `aux_scripts/filter_secretome_dtm2.py` | `05_secretome/runMe_deep.sh` |
| Effectors | EffectorP 3.0 | `05_secretome/runMe_effectorp.sh` |
| Functional annotation | InterProScan, antiSMASH, funannotate annotate | `06_functional_annotation/` |
| Genome similarity | skani | `07_genome_similarity/runMe_skani.sh` |

### References

The reference genomes go through the same QC, gene prediction, secretome and
functional annotation steps. OrthoFinder is run on the reference set
(`06_functional_annotation/runMe_funannotate_orthofinder.sh`).

### Phylogenomics

1. BUSCO (`fungi_odb12`) is run on the isolates and the reference genomes
   (`runMe_busco.sh`, or as a single job with `qsubMe_busco.sh`).
2. `runMe_busco_phylo.sh` collects the results for the retained isolates and
   the references, builds supermatrices at 70, 80 and 90% single-copy
   occupancy with BUSCO_phylogenomics, and labels the reference tips with
   their taxonomy.
3. `qsubMe_iqtree_unpartitioned.sh` infers the tree from the 70% supermatrix
   with IQ-TREE (`-m Q.YEAST+F+I+R10 -B 1000`). The tree includes a
   Basidiomycota outgroup and is rooted on it in the Quarto document.

### Parsers (`aux_scripts/`)

| Script | Input | Output |
|---|---|---|
| `filter_secretome_dtm2.py` | SignalP + DeepTMHMM2 results | TM-filtered secretome used as EffectorP input |
| `parse_effectorp.py` | EffectorP 3.0 TSV | Per-genome effector counts by class (apoplastic, cytoplasmic, dual) |
| `parse_antismash.py` | antiSMASH JSON | Per-genome BGC class counts and a long table of regions |
| `parse_funannotate_gff3_annotations.py` | funannotate GFF3 | CAZy, COG, EC and product summaries, normalised per 1000 genes |
| `parse_merqury_qv.py` | Merqury `.qv` files | Per-isolate consensus quality table |

Each script documents its input format and command-line options in its header.

## Statistical analysis and figures

`R/comparative_fungi_code.qmd` has two parts:

1. **Data processing and statistical analyses:**
   - genome QC filtering (BUSCO complete ≥ 70%, duplicated < 15%)
   - taxonomic placement of each isolate from the tree
   - PERMANOVA and betadisper on Hellinger-transformed CAZy, secreted CAZy,
     effector, COG and BGC profiles, and on a combined block-weighted matrix
   - matched comparisons of Dothideomycete endophytes against the needle
     pathogens
   - Mantel tests between phylogenetic and functional distances
2. **Figures:** every chunk whose label starts with `fig_`. These chunks reuse
   objects created in part 1, so the document must be run from top to bottom.

The document reads its input tables from `files/` (relative to the `.qmd`)
and writes figures to `figures/`. Most tables are produced by
`aux_scripts/parse_funannotate_gff3_annotations.py` and the other parsers.
The others are `busco_summary.tsv`, `quast_summary.tsv`,
`repeatmasker_summary.csv`, `references.csv` (reference metadata) and the
IQ-TREE tree file.

R packages: tidyverse, patchwork, paletteer, ggrepel, ggpubr, ggtext, ggpmisc,
scales, pheatmap, ComplexHeatmap, circlize, vegan, ape, phangorn.

## Running the scripts

The shell scripts were written for a PBS cluster. Each `runMe_*.sh` script
loops over the samples, skips those that are already finished or running,
and submits one `qsub` job per sample. The `qsubMe_*.sh` scripts are
single PBS jobs. Every tool runs in its own conda environment.

**The scripts are not portable as they stand.** They use absolute paths on
the original cluster. They also `source` a configuration file
(`sourceMe_configs.sh`) that defines the sample lists, directories and helper
functions, and that file is not in this repository. To rerun the pipeline,
adapt the paths and provide the equivalent variables (for example
`FILTERED_SAMPLES`, `LOWDUP_SAMPLES`, `HOMEDIR`, `SCRATCHDIR` and
`AUX_SCRIPTS`).

The Python parsers need Python 3 and pandas. The other modules they use are
from the standard library.

## Tool versions

| Tool | Version | Notes |
|------|---------|-------|
| Flye | 2.9.6 | |
| Medaka | 2.2.2 | Model `r1041_e82_400bps_sup_v5.2.0` |
| LINKS | 2.0.1 | |
| QUAST | 5.3.0 | |
| BUSCO | 6.1.0 | Lineage dataset `fungi_odb12` |
| Merqury | 1.4.1 | |
| RepeatMasker | 4.2.2 | |
| funannotate | 1.8.16 | |
| InterProScan | 5.77-108.0 | |
| antiSMASH | 8.0.4 | |
| SignalP | 6.0 | |
| DeepTMHMM2 | | |
| EffectorP | 3.0.0-beta | |
| skani | 0.3.2 | |
| BUSCO_phylogenomics | 20260630 | |
| IQ-TREE | 3.1.4 | |

## Data availability

- Raw reads and Genomes: BioProject PRJNA1474061

## Citation

If you use this work, please cite:

