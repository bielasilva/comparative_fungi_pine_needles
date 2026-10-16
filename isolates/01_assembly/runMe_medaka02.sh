#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/medaka02 || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/medaka_pbs_jobid.csv ]]; then
    mkdir -p logs/medaka02
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/medaka_pbs_jobid.csv
else
    cp logs/pbs_jobid/medaka_pbs_jobid.csv logs/pbs_jobid/medaka_pbs_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} ))
echo "Total checks to consider: ${total_checks}"

FINAL_DIR="${HOMEDIR}/results/assemblies/medaka/02"
mkdir -p $FINAL_DIR

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    NANOPORE_FQ="${SCRATCHDIR}/fastq/long_reads/${SAMPLE}.fastq.gz"
    MEDAKA_FASTA="${HOMEDIR}/results/assemblies/medaka/01/${SAMPLE}.medaka01.fasta"
    OUT_DIR="${SCRATCHDIR}/results/assemblies/medaka/02/${SAMPLE}"

    # Check if it has been run before
    if [[ ! -s ${NANOPORE_FQ} || ! -s ${MEDAKA_FASTA} ]]; then
        progress_log "red" "medaka for ${SAMPLE} cannot run because the input file does not exist. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}.medaka02.fasta" ]]; then
        progress_log "green" "Medaka for ${SAMPLE} already completed. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/consensus.fasta" ]]; then
        progress_log "green" "Medaka for ${SAMPLE} already completed, but not copied. Copying."
        cp ${OUT_DIR}/consensus.fasta ${FINAL_DIR}/${SAMPLE}.medaka02.fasta
        continue
    else
        if qstat -f | grep -q "Job_Name = ${SAMPLE}_medaka02"; then
            progress_log "yellow" "medaka for ${SAMPLE} is already running. Skipping."
            continue
        elif qstat -f | grep -q -e "Job_Name = ${SAMPLE}_medaka01"; then
            progress_log "red" "Medaka01 for ${SAMPLE} is still running. Skipping."
            continue
        else
            wait_checker
            progress_log "magenta" "Submitting medaka for ${SAMPLE}"
            new_jobs=$((new_jobs + 1))
            sed -i "/${SAMPLE}_medaka02/d" logs/pbs_jobid/medaka_pbs_jobid.csv
            rm -rf ${OUT_DIR}
            mkdir -p ${OUT_DIR}
        fi
    fi

    qsub <<- EOF | sed -e "s/^/${SAMPLE}_medaka02,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/medaka_pbs_jobid.csv
#!/bin/env bash

#PBS -N ${SAMPLE}_medaka02
#PBS -j oe
#PBS -o logs/medaka02/${SAMPLE}_medaka02.out
#PBS -l walltime=24:00:00
#PBS -l nodes=1
#PBS -l ncpus=20
#PBS -l mem=120gb
#PBS -q large

source /home/aubgxs001/miniforge3/bin/activate Medaka

# Run medaka
medaka_consensus \
    -t 20 \
    -i ${NANOPORE_FQ} \
    -d ${MEDAKA_FASTA} \
    -m r1041_e82_400bps_sup_v5.2.0 \
    -o ${OUT_DIR}


if [[ ! -s ${OUT_DIR}/consensus.fasta ]]; then
    echo "medaka failed for ${SAMPLE}. Exiting."
    exit 1
else
    ln -sf ${OUT_DIR}/consensus.fasta ${OUT_DIR}/${SAMPLE}_medaka02.fasta
fi

EOF

progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

done


progress_done "Done submitting medaka on all samples. Submitted ${new_jobs} new jobs."