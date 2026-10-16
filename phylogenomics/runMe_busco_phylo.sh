#!/bin/env bash

# Define general configs
config_walltime=05:00:00
config_ncpus=30
config_mem=20gb
config_queue=large

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

BUSCO_PHYLO_DIR="/home/aubgxs001/scratch/comparative_fungi/busco_phylo_links"
SUMMARY="/home/aubgxs001/scratch/databases/reference_genomes/assembly_summary_all_taxonomy_2026-09-09.tsv"
GENOMES_DIR="${BUSCO_PHYLO_DIR}/genomes"

mkdir -p ${BUSCO_PHYLO_DIR}/filtered_results

for SAMPLE in "${LOWDUP_SAMPLES[@]}"; do
    ln -s -f "${BUSCO_PHYLO_DIR}/results/${SAMPLE}" ${BUSCO_PHYLO_DIR}/filtered_results
done

REFERENCES=($(awk 'NR > 1 {print $1}' /home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/ref_list.txt | sort -u))
for REFERENCE in ${REFERENCES[@]}; do
    ln -s -f ${BUSCO_PHYLO_DIR}/results/${REFERENCE} ${BUSCO_PHYLO_DIR}/filtered_results
done

for PATHOGEN_PATH in $(ls /home/aubgxs001/scratch/comparative_fungi/references/*/*.fna);do
    PATHOGEN_FULL=$(basename $PATHOGEN_PATH)
    PATHOGEN="${PATHOGEN_FULL%%.*}"
    ln -s -f ${BUSCO_PHYLO_DIR}/results/${PATHOGEN} ${BUSCO_PHYLO_DIR}/filtered_results
done

for i in "90" "80" "70"; do

    OUT_DIR="${BUSCO_PHYLO_DIR}/busco_phylo_${i}_outgroup"
    # rm -r $OUT_DIR

bash <<- EOF
#!/bin/env bash

#PBS -N busco_phylo_${i}_outgroup
#PBS -j oe
#PBS -o logs/busco_phylo_${i}_outgroup.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

# source /home/aubgxs001/miniforge3/bin/activate BUSCO_phylo

# Run BUSCO

# BUSCO_phylogenomics.py \
#     --input ${BUSCO_PHYLO_DIR}/filtered_results \
#     --output ${OUT_DIR} \
#     --threads 30 \
#     --percent_single_copy ${i}


echo "Doing $i"
cp ${OUT_DIR}/gene_trees/ALL.tree ${OUT_DIR}/gene_trees/phylo_${i}_tax.tree

references=(\$(awk 'NR > 1 {print \$1}' /home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/ref_list.txt | sort -u))
for reference in \${references[@]}; do
    reference_stpd="\${reference%%.*}"
    taxonomy=\$(awk -F'\t' -v OFS='|' -v ref="\$reference_stpd" '\$1 ~ ref {print \$40, \$41, \$42, \$43, \$45}' "$SUMMARY")
    sed -i "s/\$reference_stpd/\$taxonomy/g" ${OUT_DIR}/gene_trees/phylo_${i}_tax.tree
done
EOF
done