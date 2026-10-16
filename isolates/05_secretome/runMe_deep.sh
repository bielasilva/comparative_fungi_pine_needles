#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=24:00:00
config_ncpus=5
config_mem=15gb
config_queue=gpu

# Create output directories
if [[ ! -d logs/deeptmhmm || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/deeptmhmm_pbs_jobid.csv ]]; then
    mkdir -p logs/deeptmhmm
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/deeptmhmm_pbs_jobid.csv
else
    cp logs/pbs_jobid/deeptmhmm_pbs_jobid.csv logs/pbs_jobid/deeptmhmm_pbs_jobid.csv.bk
fi

total_checks=$(( ${#FILTERED_SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${FILTERED_SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

    FUN_DIR="${SCRATCHDIR}/results/annotation/funannotate_links/${SAMPLE}"
    SIGNALP_FA="${FUN_DIR}/signalp/processed_entries.fasta"
    OUT_DIR="${FUN_DIR}/deeptmhmm"

    # Sanity checks
    if [[ ! -s "${SIGNALP_FA}" ]]; then
        progress_log "red" "deeptmhmm for ${SAMPLE} skipped - no SignalP file found"
        continue
    elif [[ -s "${OUT_DIR}/filtered_signalp.fasta" ]]; then
        progress_log "green" "deeptmhmm for ${SAMPLE} has already completed. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = deeptmhmm_${SAMPLE}_links"; then
        progress_log "yellow" "deeptmhmm for ${SAMPLE} is already running. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = ${SAMPLE}_signalp_links"; then
        progress_log "yellow" "SignalP for ${SAMPLE} is still running. Skipping."
        continue
    else
        wait_checker
        progress_log "magenta" "Submitting deeptmhmm for ${SAMPLE}"
        new_jobs=$((new_jobs + 1))
        sed -i "/deeptmhmm_${SAMPLE}_links,/d" logs/pbs_jobid/deeptmhmm_pbs_jobid.csv
        rm -r "${OUT_DIR}" 2> /dev/null
    fi

    qsub > >(sed -e "s/^/deeptmhmm_${SAMPLE}_links,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/deeptmhmm_pbs_jobid.csv) 2> >(raise_error) <<- EOF 

#PBS -N deeptmhmm_${SAMPLE}_links
#PBS -j oe
#PBS -o logs/deeptmhmm/deeptmhmm_${SAMPLE}_links.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Dtm2

dtm2 ${SIGNALP_FA} ${OUT_DIR} --device cuda

if [[ -s ${OUT_DIR}/predicted_topologies.3line ]]; then
    python3 ${AUX_SCRIPTS}/filter_secretome_dtm2.py \
        --predictions ${OUT_DIR}/predictions.json \
        --fasta ${SIGNALP_FA} \
        --mature-input \
        --out-fasta  ${OUT_DIR}/filtered_signalp.fasta \
        --out-table  ${OUT_DIR}/secretome_TM_filter_summary.tsv
else
    echo "Something went wrong for deeptmhmm_${SAMPLE}_links. Check it"
    $(cleanup)
    exit 1
fi

$(cleanup)

EOF

done

echo_overwrite_2 "Done submitting deeptmhmm on all samples. Submitted ${new_jobs} new jobs."