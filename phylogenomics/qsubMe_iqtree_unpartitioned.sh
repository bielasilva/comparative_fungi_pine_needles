#!/bin/env bash

#PBS -N iqtree_outgroup_unpart_70
#PBS -j oe
#PBS -o logs/iqtree_outgroup_unpart_70.out
#PBS -l walltime=320:00:00
#PBS -l nodes=1
#PBS -l ncpus=30
#PBS -l mem=120gb
#PBS -q large

source /home/aubgxs001/miniforge3/bin/activate Iqtree

# Run IQ-Tree
BUSCO_PHYLO_DIR="/home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/busco_phylo_70_outgroup"
SUPERMATRIX="${BUSCO_PHYLO_DIR}/supermatrix/SUPERMATRIX.phylip"
PARTITIONS="${BUSCO_PHYLO_DIR}/supermatrix/SUPERMATRIX.partitions.nex"

mkdir -p ${BUSCO_PHYLO_DIR}/iqtree

iqtree -s $SUPERMATRIX \
    -m Q.YEAST+F+I+R10 \
    -B 1000 \
    -T AUTO -ntmax 30 \
    --prefix ${BUSCO_PHYLO_DIR}/iqtree/iqtree_outgroup_unpart_70