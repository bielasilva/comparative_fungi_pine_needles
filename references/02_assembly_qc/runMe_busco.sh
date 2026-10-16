#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=05:00:00
config_ncpus=10
config_mem=15gb
config_queue=medium

# Create output directories
if [[ ! -d logs/busco || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/busco_pbs_jobid.csv ]]; then
    mkdir -p logs/busco
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/busco_pbs_jobid.csv
else
    cp logs/pbs_jobid/busco_pbs_jobid.csv logs/pbs_jobid/busco_pbs_jobid.csv.bk
fi

FINAL_DIR="${HOMEDIR}/results/busco/references"
mkdir -p $FINAL_DIR

total_checks=$((${#REFERENCES_FULL[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE_FULL in "${REFERENCES_FULL[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"


    ASSEMBLY_FASTA="${SCRATCHDIR}/references/${SAMPLE_FULL}/${SAMPLE_FULL}_genomic.fna"
    SAMPLE=$(echo "${SAMPLE_FULL}" | tr '_' '-' | awk -F'-' '{print substr($1,1,1) $2}')
    OUT_DIR="${SCRATCHDIR}/results/busco/references"

    # Check if it has been run before
    if [[ ! -s ${ASSEMBLY_FASTA} ]]; then
        progress_log "red" "busco for ${SAMPLE} cannot run because the input files does not exist. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}.busco.json" ]]; then
        progress_log "green" "busco for ${SAMPLE} already completed. Skipping"
        # rm -r ${OUT_DIR}/${SAMPLE} 2> /dev/null
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json" ]]; then
        progress_log "green" "busco for ${SAMPLE} already completed, but not copied. Copying."
        cp ${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json ${FINAL_DIR}/${SAMPLE}.busco.json
        continue
    else
        if qstat -f | grep -q -E "Job_Name = busco_${SAMPLE}$"; then
            progress_log "yellow" "busco for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            progress_log "magenta" "Submitting busco Fungi for ${SAMPLE}."
            new_jobs=$((new_jobs + 1))
            sed -i "/busco_${SAMPLE},/d" logs/pbs_jobid/busco_pbs_jobid.csv
            rm -rf ${OUT_DIR}/${SAMPLE}
        fi
    fi

    qsub <<- EOF | sed -e "s/^/busco_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/busco_pbs_jobid.csv
#!/bin/env bash

#PBS -N busco_${SAMPLE}
#PBS -j oe
#PBS -o logs/busco/busco_${SAMPLE}.out
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
    echo "busco for ${SAMPLE} failed. Please check the log file for details: logs/busco/busco_${SAMPLE}.out"
    exit 1
else
    echo "BUSCO successful, copying data."
    cp ${OUT_DIR}/${SAMPLE}/short_summary.specific.fungi_odb12.${SAMPLE}.json ${FINAL_DIR}/${SAMPLE}.busco.json

    if [[ -s ${FINAL_DIR}/${SAMPLE}.busco.json ]]; then
        echo "Cleaning up scratch"
        rm -r ${OUT_DIR}/${SAMPLE}
    else
        echo "Copying failed, check manually"
    fi
fi
EOF

done


progress_done "Done submitting Busco on all samples. Submitted ${new_jobs} new jobs."