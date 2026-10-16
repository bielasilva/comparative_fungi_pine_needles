#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# General Configs
config_walltime="1:00:00"
config_nodes="1"
config_ncpus="4"
config_mem="4gb"
config_queue="small"

# Create output directories
if [[ ! -d logs/quast || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/quast_pbs_jobid.csv ]]; then
    mkdir -p logs/quast
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/quast_pbs_jobid.csv
else
    cp logs/pbs_jobid/quast_pbs_jobid.csv logs/pbs_jobid/quast_pbs_jobid.csv.bk
fi

FINAL_DIR="${HOMEDIR}/results/quast/references"
mkdir -p $FINAL_DIR

total_checks=$((${#REFERENCES_FULL[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE_FULL in "${REFERENCES_FULL[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"


    GENOME="${SCRATCHDIR}/references/${SAMPLE_FULL}/${SAMPLE_FULL}_genomic.fna"
    SAMPLE=$(echo "${SAMPLE_FULL}" | tr '_' '-' | awk -F'-' '{print substr($1,1,1) $2}')
    OUT_DIR="${SCRATCHDIR}/results/quast/references/${SAMPLE}"

    # Check if it has been run before
    if [[ ! -s ${GENOME} ]]; then
        progress_log "red" "QUAST for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}_report.tsv" ]]; then
        progress_log "green" "QUAST for ${SAMPLE} already completed. Skipping."
        # rm -r ${OUT_DIR}
        continue
    elif [[ -s "${OUT_DIR}/report.tsv" ]]; then
        progress_log "green" "QUAST for ${SAMPLE} already completed, but not copied. Copying."
        continue
    else
        if qstat -f | grep -q -E "Job_Name = quast_${SAMPLE}$"; then
            progress_log "yellow" "QUAST for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            progress_log "magenta" "Submitting QUAST for ${SAMPLE}."
            new_jobs=$((new_jobs + 1))
            sed -i "/quast_${SAMPLE},/d" logs/pbs_jobid/quast_pbs_jobid.csv
            # rm -rf ${OUT_DIR}
            mkdir -p ${OUT_DIR}
        fi
    fi

    # qsub <<- EOF | sed -e "s/^/quast_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/quast_pbs_jobid.csv
    qsub > >(sed -e "s/^/quast_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/quast_pbs_jobid.csv) 2> >(raise_error) <<- EOF
#!/bin/env bash

#PBS -N quast_${SAMPLE}
#PBS -j oe
#PBS -o logs/quast/quast_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

source /home/aubgxs001/miniforge3/bin/activate Quast

$(log_asax)

# Run Quast

quast.py -o ${OUT_DIR} \
    -t 4 \
    --circos \
    --min-contig 500 \
    --fungus \
    --k-mer-stats \
    ${GENOME}

if [[ ! -s "${OUT_DIR}/report.tsv" ]]; then
    echo "QUAST for ${SAMPLE} at ${PROGRAM} failed. Check the output log for details."
    exit 1
else
    cp "${OUT_DIR}/report.tsv"            "${FINAL_DIR}/${SAMPLE}_report.tsv"
    cp "${OUT_DIR}/transposed_report.tsv" "${FINAL_DIR}/${SAMPLE}_transposed_report.tsv"
    if [[ -s "${FINAL_DIR}/${SAMPLE}_report.tsv" ]]; then
        rm -r ${OUT_DIR}
    fi
fi

$(cleanup)
EOF
done


progress_done "Done submitting QUAST on all samples. Submitted ${new_jobs} new jobs."