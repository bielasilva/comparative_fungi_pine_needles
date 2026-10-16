#!/bin/env bash

#PBS -N busco_run
#PBS -j oe
#PBS -o logs/busco_run.out
#PBS -l walltime=24:00:00
#PBS -l nodes=1
#PBS -l ncpus=50
#PBS -l mem=50gb
#PBS -q large

source /home/aubgxs001/miniforge3/bin/activate Busco

# Run BUSCO

GENOMES_DIR="/home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/genomes"

mkdir -p $GENOMES_DIR

OUT_DIR="/home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/results"

# rm ${GENOMES_DIR}/assembly_summary.csv 2> /dev/null

# ln -s -f /home/aubgxs001/comparative_fungi/results/assemblies/repeatmasker/links/*.masked.fasta $GENOMES_DIR

# references=($(awk 'NR > 1 {print $1}' /home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/ref_list.txt | sort -u))
# for reference in ${references[@]}; do
#     ln -s -f /home/aubgxs001/scratch/databases/reference_genomes/fasta/${reference}_genomic.fna.gz ${GENOMES_DIR}/
#     grep $reference /home/aubgxs001/scratch/databases/reference_genomes/assembly_summary_contig_no-filtering_subset-genus.csv >> ${GENOMES_DIR}/assembly_summary.csv
# done

# for genome in $(ls $GENOMES_DIR); do
#     if [[ -s "/home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/BUSCO_genomes/$genome/run_fungi_odb12/short_summary.json" ]]; then
#         echo "$genome already done"
#     else
#         busco --in ${GENOMES_DIR}/$genome -f \
#             --lineage_dataset /scratch/aubgxs001/databases/busco_downloads/lineages/fungi_odb12 \
#             --cpu 20 \
#             --mode genome \
#             --out_path ${OUT_DIR}/BUSCO_genomes/
#     fi
# done

busco --in ${GENOMES_DIR} -f \
    --lineage_dataset /scratch/aubgxs001/databases/busco_downloads/lineages/fungi_odb12 \
    --cpu 50 \
    --mode genome \
    --out_path ${OUT_DIR}