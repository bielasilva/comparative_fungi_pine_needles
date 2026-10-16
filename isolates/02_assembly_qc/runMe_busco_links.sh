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

total_checks=$(( ${#SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

FINAL_DIR="${HOMEDIR}/results/busco/links"
mkdir -p $FINAL_DIR

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

    ASSEMBLY_FASTA="${HOMEDIR}/results/assemblies/links/${SAMPLE}.links.fasta"
    OUT_DIR="${SCRATCHDIR}/results/busco/links"

    # Check if it has been run before
    if [[ ! -s ${ASSEMBLY_FASTA} ]]; then
        progress_log "red" "busco_links for ${SAMPLE} cannot run because the input files does not exist. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}.busco_links.json" ]]; then
        progress_log "green" "busco_links for ${SAMPLE} already completed. Skipping"
        # rm -r ${OUT_DIR}/${SAMPLE} 2> /dev/null
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json" ]]; then
        progress_log "green" "busco_links for ${SAMPLE} already completed, but not copied. Copying."
        cp ${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json ${FINAL_DIR}/${SAMPLE}.busco_links.json
        continue
    else
        if qstat -f | grep -q -E "Job_Name = busco_links_${SAMPLE}$"; then
            progress_log "yellow" "busco_links for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            new_jobs=$((new_jobs + 1))
            progress_log "magenta" "Submitting busco_links Fungi for ${SAMPLE}."
            sed -i "/busco_links_${SAMPLE},/d" logs/pbs_jobid/busco_links_pbs_jobid.csv
            rm -rf ${OUT_DIR}/${SAMPLE}
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

# Run BUSCO

busco --in ${ASSEMBLY_FASTA} -f \
    --out ${SAMPLE} \
    --lineage_dataset /scratch/aubgxs001/databases/busco_downloads/lineages/fungi_odb12  \
    --cpu 10 \
    --mode genome \
    --out_path ${OUT_DIR}

if [[ ! -s ${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json ]]; then
    echo "busco_links for ${SAMPLE} failed. Please check the log file for details: logs/busco_links/busco_links_${SAMPLE}.out"
    exit 1
else
    echo "BUSCO successful, copying data."
    cp ${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json ${FINAL_DIR}/${SAMPLE}.busco_links.json

    if [[ -s ${FINAL_DIR}/${SAMPLE}.busco_links.json ]]; then
        echo "Cleaning up scratch"
        rm -r ${OUT_DIR}/${SAMPLE}
    else
        echo "Copying failed, check manually"
    fi
fi
EOF

done


progress_done "Done submitting Busco on all samples. Submitted ${new_jobs} new jobs."