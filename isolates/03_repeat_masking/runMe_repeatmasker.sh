#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Create output directories
if [[ ! -d logs/repeatmasker || ! -d logs/pbs_jobid || -s logs/pbs_jobid/repeatmasker_pbs_jobid.csv ]]; then
    mkdir -p logs/repeatmasker
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/repeatmasker_pbs_jobid.csv
else
    cp logs/pbs_jobid/repeatmasker_pbs_jobid.csv logs/pbs_jobid/repeatmasker_pbs_jobid.csv.bk
fi

total_checks=$(( ${#SHORTREAD_SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

FINAL_DIR="${HOMEDIR}/results/assemblies/repeatmasker/polypolish"
mkdir -p $FINAL_DIR

for SAMPLE in "${SHORTREAD_SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    # FASTA="${SCRATCHDIR}/fastas/${SAMPLE}.fasta"
    FASTA="${HOMEDIR}/results/assemblies/polypolish/${SAMPLE}.polypolish.fasta"
    OUT_DIR="${SCRATCHDIR}/results/repeatmasker_polypolish/${SAMPLE}"

    # Check if it has been run before
    if [[ ! -s ${FASTA} ]]; then
        progress_log "red" "Repeatmasker for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}.masked.fasta" ]]; then
        progress_log "green" "polypolish for ${SAMPLE} already completed. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}.cut.masked" ]]; then
        progress_log "green" "Repeatmasker for ${SAMPLE} already completed but not copied. Copying."
        cp ${OUT_DIR}/${SAMPLE}.cut.masked ${FINAL_DIR}/${SAMPLE}.masked.fasta
        rm -r ${OUT_DIR}
        continue
    else
        if qstat -f | grep -q -E "Job_Name = repeatmasker_${SAMPLE}_polypolish"; then
            progress_log "yellow" "repeatmasker for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            new_jobs=$((new_jobs + 1))
            progress_log "magenta" "Submitting repeatmasker for ${SAMPLE}"
            sed -i "repeatmasker_${SAMPLE}_polypolish,/d" logs/pbs_jobid/repeatmasker_pbs_jobid.csv
            rm -rf ${OUT_DIR}
            mkdir -p ${OUT_DIR}
        fi
    fi

    # cat <<- EOF 
    qsub <<- EOF | sed -e "s/^/repeatmasker_${SAMPLE}_polypolish,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/repeatmasker_pbs_jobid.csv
#!/bin/env bash

#PBS -N repeatmasker_${SAMPLE}_polypolish
#PBS -j oe
#PBS -o logs/repeatmasker/repeatmasker_${SAMPLE}_polypolish.out
#PBS -l walltime=2:00:00
#PBS -l nodes=1
#PBS -l ncpus=15
#PBS -l mem=15gb
#PBS -q medium

source /home/aubgxs001/miniforge3/bin/activate Repeatmasker

# Run repeatmasker

cut -d ',' -f1 ${FASTA} > ${OUT_DIR}/${SAMPLE}.cut

RepeatMasker -pa 15 -species fungi -xsmall -dir ${OUT_DIR} ${OUT_DIR}/${SAMPLE}.cut

if [[ ! -s ${OUT_DIR}/${SAMPLE}.cut.masked ]]; then
    echo "Repeatmasker failed for ${SAMPLE}. Exiting."
    exit 1
else
    echo "Repeatmasker successful, copying to home."
    cp ${OUT_DIR}/${SAMPLE}.cut.masked ${FINAL_DIR}/${SAMPLE}.masked.fasta
    rm -r ${OUT_DIR}
fi

EOF
done

echo_overwrite_2 "Done submitting repeatmasker on all samples. Submitted ${new_jobs} new jobs."