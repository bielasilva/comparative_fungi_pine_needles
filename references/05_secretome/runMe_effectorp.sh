#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=01:00:00
config_ncpus=4
config_mem=4gb
config_queue=small

# Create output directories
if [[ ! -d logs/effectorp || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/effectorp_pbs_jobid.csv ]]; then
    mkdir -p logs/effectorp
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/effectorp_pbs_jobid.csv
else
    cp logs/pbs_jobid/effectorp_pbs_jobid.csv logs/pbs_jobid/effectorp_pbs_jobid.csv.bk
fi

# Read the list of samples from the CSV file, skipping the header and make spaces in the name underscores
readarray -t REFERENCES < <(awk -F',' 'NR > 1 {print $1}' /home/aubgxs001/scratch/comparative_fungi/references/reference_genomes.csv)
REFERENCES=("${REFERENCES[@]// /_}")

total_checks=$(( ${#REFERENCES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE_FULL in "${REFERENCES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    SAMPLE=$(echo "${SAMPLE_FULL}" | tr '_' '-' | awk -F'-' '{print substr($1,1,1) $2}')

    FUN_DIR="${SCRATCHDIR}/results/annotation/funannotate/${SAMPLE}"
    SIGNALP_FA="${FUN_DIR}/deeptmhmm/${SAMPLE}_filtered_signalp.fasta"
    OUT_DIR="${FUN_DIR}/effectorp"

    # Sanity checks
    if [[ ! -s "${SIGNALP_FA}" ]]; then
        progress_log "red" "EffectorP for ${SAMPLE} skipped - no SignalP file found"
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}.tsv" ]]; then
        progress_log "green" "EffectorP for ${SAMPLE} has already completed. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = effectorp_${SAMPLE}"; then
        progress_log "yellow" "EffectorP for ${SAMPLE} is already running. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = ${SAMPLE}_signalp"; then
        progress_log "yellow" "SignalP for ${SAMPLE} is still running. Skipping."
        continue
    else
        wait_checker
        progress_log "magenta" "Submitting EffectorP for ${SAMPLE}"
        new_jobs=$((new_jobs + 1))
        sed -i "/effectorp_${SAMPLE},/d" logs/pbs_jobid/effectorp_pbs_jobid.csv
        mkdir -p "${OUT_DIR}"
    fi

    qsub > >(sed -e "s/^/effectorp_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/effectorp_pbs_jobid.csv) 2> >(raise_error) <<- EOF

#PBS -N effectorp_${SAMPLE}
#PBS -j oe
#PBS -o logs/effectorp/effectorp_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

cd "${OUT_DIR}"

EffectorP="python3 /home/aubgxs001/scratch/tools/EffectorP_3.0.0-beta/EffectorP.py"

Seqkit="mamba run -n Seqkit seqkit"

\$EffectorP -f -i ${SIGNALP_FA} -o "${OUT_DIR}/${SAMPLE}.tsv"

# retained IDs from the filter output
cut -f1 ${FUN_DIR}/deeptmhmm/secretome_TM_filter_summary.tsv | tail -n +2 > ${FUN_DIR}/deeptmhmm/kept_ids.txt

# pull full-length versions from the original proteome
\$Seqkit grep -f ${FUN_DIR}/deeptmhmm/kept_ids.txt ${SIGNALP_FA} > ${OUT_DIR}/${SAMPLE}_secretome_fulllength.faa

\$EffectorP -f -i ${OUT_DIR}/${SAMPLE}_secretome_fulllength.faa -o ${OUT_DIR}/${SAMPLE}_full.tsv

$(cleanup)
EOF

done

echo_overwrite_2 "Done submitting EffectorP on all samples. Submitted ${new_jobs} new jobs."