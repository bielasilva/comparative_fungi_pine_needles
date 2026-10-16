#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=06:00:00
config_ncpus=8
config_mem=15gb
config_queue=medium

# Create output directories
if [[ ! -d logs/merqury || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/merqury_pbs_jobid.csv ]]; then
    mkdir -p logs/merqury
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/merqury_pbs_jobid.csv
else
    cp logs/pbs_jobid/merqury_pbs_jobid.csv logs/pbs_jobid/merqury_pbs_jobid.csv.bk
fi

total_checks=$(( ${#SHORTREAD_SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${SHORTREAD_SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

    READS_R1="${SCRATCHDIR}/fastq/short_reads/fastp/${SAMPLE}_trimmed_R1.fastq.gz"
    READS_R2="${SCRATCHDIR}/fastq/short_reads/fastp/${SAMPLE}_trimmed_R2.fastq.gz"
    ASM_FASTA="${SCRATCHDIR}/fastas/${SAMPLE}.fasta"   # EDIT: confirm filename pattern

    MERQ_DIR="${SCRATCHDIR}/results/merqury/${SAMPLE}"
    OUT_DIR="${MERQ_DIR}/work"
    FINAL_DIR="${HOMEDIR}/results/merqury"
    MERQURY_K=31

    # Sanity checks
    if [[ ! -s "${READS_R1}" || ! -s "${READS_R2}" ]]; then
        progress_log "yellow" "Merqury for ${SAMPLE} skipped - no matched Illumina reads found."
        continue
    elif [[ ! -s "${ASM_FASTA}" ]]; then
        progress_log "red" "Merqury for ${SAMPLE} cannot be run because the assembly is not available. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}.qv" ]]; then
        progress_log "green" "Merqury for ${SAMPLE} has already completed and moved. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = merqury_${SAMPLE}"; then
        progress_log "yellow" "Merqury for ${SAMPLE} is already running. Skipping."
        continue
    else
        wait_checker
        new_jobs=$((new_jobs + 1))
        progress_log "magenta" "Submitting Merqury for ${SAMPLE}"
        sed -i "/merqury_${SAMPLE},/d" logs/pbs_jobid/merqury_pbs_jobid.csv
        mkdir -p "${OUT_DIR}"
    fi

    # qsub <<- EOF 2> >(progress_log_error) | sed -e "s/^/merqury_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/merqury_pbs_jobid.csv
    qsub 2> >(raise_error) <<- EOF > >(sed -e "s/^/merqury_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/merqury_pbs_jobid.csv)

#PBS -N merqury_${SAMPLE}
#PBS -j oe
#PBS -o logs/merqury/merqury_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Merqury

cd "${OUT_DIR}"

# Build meryl k-mer dbs per mate, then merge - counting R1 and R2 in one
# 'meryl count' call does NOT combine both mates' counts, they must be
# built separately and merged with union-sum (marbl/merqury wiki, "Prepare
# meryl dbs"). Skipped if a previous partial run already left the merged db.
if [[ ! -d "${SAMPLE}.reads.meryl" ]]; then
    meryl k=${MERQURY_K} count threads=${config_ncpus} output ${SAMPLE}.R1.meryl "${READS_R1}"
    meryl k=${MERQURY_K} count threads=${config_ncpus} output ${SAMPLE}.R2.meryl "${READS_R2}"
    meryl union-sum output ${SAMPLE}.reads.meryl ${SAMPLE}.R1.meryl ${SAMPLE}.R2.meryl
    # rm -rf ${SAMPLE}.R1.meryl ${SAMPLE}.R2.meryl
fi

merqury.sh ${SAMPLE}.reads.meryl "${ASM_FASTA}" ${SAMPLE}

mkdir -p ${FINAL_DIR}
cp ${OUT_DIR}/${SAMPLE}.qv ${OUT_DIR}/${SAMPLE}.completeness.stats ${FINAL_DIR}/ 2>/dev/null
cp ${OUT_DIR}/${SAMPLE}*.spectra-cn.*.png ${FINAL_DIR}/ 2>/dev/null

if [[ -s ${FINAL_DIR}/${SAMPLE}.qv ]]; then
    echo "Merqury for ${SAMPLE} completed successfully."
    # rm -rf ${MERQ_DIR}/work ${SAMPLE}.*.meryl ${SAMPLE}.${SAMPLE}.meryl 2>/dev/null
else
    echo "Merqury for ${SAMPLE} did not complete successfully. Please check the logs for more details."
fi
EOF

done

echo_overwrite_2 "Done submitting Merqury on all samples. Submitted ${new_jobs} new jobs."