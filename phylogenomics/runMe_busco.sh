#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=05:00:00
config_ncpus=10
config_mem=15gb
config_queue=medium

# Create output directories
if [[ ! -d logs/busco_links || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/busco_links_pbs_jobid.csv ]]; then
    mkdir -p logs/busco_links
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/busco_links_pbs_jobid.csv
else
    cp logs/pbs_jobid/busco_links_pbs_jobid.csv logs/pbs_jobid/busco_links_pbs_jobid.csv.bk
fi

GENOMES_DIR="${SCRATCHDIR}/busco_phylo_links/genomes"
mkdir -p $GENOMES_DIR

OUT_DIR="${SCRATCHDIR}/busco_phylo_links/results"

# ln -s -f /home/aubgxs001/comparative_fungi/results/assemblies/repeatmasker/links/*.masked.fasta $GENOMES_DIR
# ln -s -f /home/aubgxs001/scratch/comparative_fungi/references/*/*.fna $GENOMES_DIR
references=($(awk 'NR > 1 {print $1}' /home/aubgxs001/scratch/comparative_fungi/busco_phylo_links/ref_list.txt | sort -u))
for reference in ${references[@]}; do
    ln -s -f /home/aubgxs001/scratch/databases/reference_genomes/fasta/${reference}*_genomic.fna.gz ${GENOMES_DIR}/
done

for SAMPLE_PATH in $(ls $GENOMES_DIR/*.f*a*); do
    
    SAMPLE_FULL=$(basename $SAMPLE_PATH)
    ASSEMBLY_FASTA="${GENOMES_DIR}/${SAMPLE_FULL}"

    SAMPLE="${SAMPLE_FULL%%.*}"

    # Check if it has been run before
    if [[ ! -s ${ASSEMBLY_FASTA} ]]; then
        echo_log "red" "busco_links for ${SAMPLE} cannot run because the input files does not exist. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json" ]]; then
        echo_log "green" "busco_links for ${SAMPLE} already completed."
        continue
    else
        if qstat -f | grep -q -E "Job_Name = busco_links_${SAMPLE}$"; then
            echo_log "yellow" "busco_links for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            echo_log "magenta" "Submitting busco_links Fungi for ${SAMPLE}."
            sed -i "/busco_links_${SAMPLE},/d" logs/pbs_jobid/busco_links_pbs_jobid.csv
            # rm -rf ${OUT_DIR}/${SAMPLE}
        fi
    fi

    qsub <<- EOF | sed -e "s/^/busco_links_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/busco_links_pbs_jobid.csv
#!/bin/env bash

#PBS -N busco_links_${SAMPLE}
#PBS -j oe
#PBS -o logs/busco_links/busco_links_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

source /home/aubgxs001/miniforge3/bin/activate Busco

$(log_asax)

# Run BUSCO

busco --in ${ASSEMBLY_FASTA} -f \
    --out ${SAMPLE} \
    --lineage_dataset /scratch/aubgxs001/databases/busco_downloads/lineages/fungi_odb12  \
    --cpu 10 \
    --mode genome \
    --out_path ${OUT_DIR}

$(cleanup)
EOF

done


progress_done "Done submitting Busco on all samples. Submitted ${new_jobs} new jobs."