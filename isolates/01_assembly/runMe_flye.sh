#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/flye || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/flye_pbs_jobid.csv ]]; then
    mkdir -p logs/flye
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/flye_pbs_jobid.csv
else
    cp logs/pbs_jobid/flye_pbs_jobid.csv logs/pbs_jobid/flye_pbs_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} ))
echo "Total checks to consider: ${total_checks}"

FINAL_DIR="${HOMEDIR}/results/assemblies/flye"
mkdir -p $FINAL_DIR

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    NANOPORE_FQ="${SCRATCHDIR}/fastq/long_reads/${SAMPLE}.fastq.gz"
    OUT_DIR="${SCRATCHDIR}/results/assemblies/flye/${SAMPLE}"

    # Check if it has been run before
    if [[ ! -s ${NANOPORE_FQ} ]]; then
        progress_log "red" "FLYE for ${SAMPLE} cannot run because the input file does not exist. Skipping."
        continue
    elif grep -q "FLYE failed for" "logs/flye/${SAMPLE}_flye.out" 2> /dev/null; then
        progress_log "red" "FLYE for ${SAMPLE} failed with error. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}.flye.fasta" ]]; then
        progress_log "green" "Flye for ${SAMPLE} already completed. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/assembly.fasta" ]]; then
        progress_log "green" "Flye for ${SAMPLE} already completed, but not copied. Copying."
        cp ${OUT_DIR}/assembly.fasta ${FINAL_DIR}/${SAMPLE}.flye.fasta
        continue
    else
        if qstat -f | grep -q "Job_Name = ${SAMPLE}_flye"; then
            progress_log "yellow" "FLYE for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            progress_log "magenta" "Submitting FLYE for ${SAMPLE}"
            new_jobs=$((new_jobs + 1))
            sed -i "/${SAMPLE}_flye/d" logs/pbs_jobid/flye_pbs_jobid.csv
            rm -rf ${OUT_DIR}
            mkdir -p ${OUT_DIR}
        fi
    fi

    qsub <<- EOF | sed -e "s/^/${SAMPLE}_flye,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/flye_pbs_jobid.csv
#!/bin/env bash

#PBS -N ${SAMPLE}_flye
#PBS -j oe
#PBS -o logs/flye/${SAMPLE}_flye.out
#PBS -l walltime=120:00:00
#PBS -l nodes=1
#PBS -l ncpus=10
#PBS -l mem=100gb
#PBS -q large

source /home/aubgxs001/miniforge3/bin/activate Flye

# Run FLYE
flye \
    --threads 10 \
    --nano-raw ${NANOPORE_FQ} \
    --min-overlap 1000 \
    --out-dir ${OUT_DIR}

if [[ ! -s ${OUT_DIR}/assembly.fasta ]]; then
    echo "FLYE failed for ${SAMPLE}. Exiting."
    exit 1
else
    echo "Flye successful, copying to home."
    cp ${OUT_DIR}/assembly.fasta ${FINAL_DIR}/${SAMPLE}.flye.fasta
fi

EOF

progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

done


progress_done "Done submitting FLYE on all samples. Submitted ${new_jobs} new jobs."