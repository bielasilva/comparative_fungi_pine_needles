#!/bin/env bash
set -e

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=12:00:00
config_ncpus=8
config_mem=10gb
config_queue=medium

# Create output directories
if [[ ! -d logs/signalp || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/signalp_pbs_jobid.csv ]]; then
    mkdir -p logs/signalp
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/signalp_pbs_jobid.csv
else
    cp logs/pbs_jobid/signalp_pbs_jobid.csv logs/pbs_jobid/signalp_pbs_jobid.csv.bk
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
    PROTEIN_FASTA="${FUN_DIR}/predict/predict_results/fungi_${SAMPLE}.proteins.fa"
    OUT_DIR="${FUN_DIR}/signalp"

    if [[ ! -s "${PROTEIN_FASTA}" ]]; then
        progress_log "red" "signalp for ${SAMPLE} cannot run because the protein file does not exist. Skipping."
        continue
    elif [[ -s "${FUN_DIR}/signalp/prediction_results.txt" ]]; then
        progress_log "green" "signalp for ${SAMPLE} already completed."
        continue
    elif [[ -s "${FUN_DIR}/signalp/${SAMPLE}/prediction_results.txt" ]]; then
        progress_log "green" "signalp for ${SAMPLE} already completed but in the wrong folder. Moving it."
        mv ${FUN_DIR}/signalp/${SAMPLE}/* ${FUN_DIR}/signalp/
        rmdir ${FUN_DIR}/signalp/${SAMPLE}
        continue
    else
        if qstat -f | grep -q "Job_Name = ${SAMPLE}_signalp"; then
            progress_log "yellow" "signalp for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            progress_log "magenta" "Submitting signalp for ${SAMPLE}"
            new_jobs=$((new_jobs + 1))
            sed -i "/${SAMPLE}_signalp/d" logs/pbs_jobid/signalp_pbs_jobid.csv
            # rm -rf ${OUT_DIR}
            mkdir -p ${OUT_DIR}
        fi
    fi
    qsub 2> >(raise_error) > >(sed -e "s/^/${SAMPLE}_signalp,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/signalp_pbs_jobid.csv) <<- EOF
#!/bin/env bash

#PBS -N ${SAMPLE}_signalp_links
#PBS -j oe
#PBS -o logs/signalp/${SAMPLE}_signalp.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

source /home/aubgxs001/miniforge3/bin/activate Signalp

$(log_asax)

cd ${OUT_DIR}

signalp6 --fastafile ${PROTEIN_FASTA} \
    --organism eukarya \
    --output_dir ${OUT_DIR} \
    --format none --mode fast

$(cleanup)
EOF

progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

done


progress_done "Done submitting signalp on all samples. Submitted ${new_jobs} new jobs."