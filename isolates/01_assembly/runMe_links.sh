#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# General Configs
config_walltime="5:00:00"
config_nodes="1"
config_ncpus="4"
config_mem="50gb"
config_queue="large"

# Create output directories
if [[ ! -d logs/links || ! -d logs/pbs_jobid || ! -s logs/pbs_jobid/links_pbs_jobid.csv ]]; then
    mkdir -p logs/links
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/links_pbs_jobid.csv
else
    cp logs/pbs_jobid/links_pbs_jobid.csv logs/pbs_jobid/links_pbs_jobid.csv.bk
fi

total_checks=$(( ${#SAMPLES[@]} ))
echo "Total checks to consider: ${total_checks}"

FINAL_DIR="${HOMEDIR}/results/assemblies/links"
mkdir -p $FINAL_DIR

for SAMPLE in "${SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    NANOPORE_FQ="${SCRATCHDIR}/fastq/long_reads/${SAMPLE}.fastq.gz"
    MEDAKA_FASTA="${SCRATCHDIR}/results/assemblies/medaka/02/${SAMPLE}/${SAMPLE}_medaka02.fasta"
    OUT_DIR="${SCRATCHDIR}/results/assemblies/links_multi/${SAMPLE}"

    # Check if it has been run before
    if [[ ! -s ${NANOPORE_FQ} || ! -s ${MEDAKA_FASTA} ]]; then
        progress_log "red" "Links for ${SAMPLE} cannot run because the input file does not exist. Skipping."
        continue
    elif [[ -s "${FINAL_DIR}/${SAMPLE}.links.fasta" ]]; then
        progress_log "green" "Links for ${SAMPLE} already completed. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}.links.scaffolds.fa" ]]; then
        progress_log "green" "Links for ${SAMPLE} already completed, but not copied. Copying."
        cp ${OUT_DIR}/${SAMPLE}.links.scaffolds.fa ${FINAL_DIR}/${SAMPLE}.links.fasta
        continue
    else
        if qstat -f | grep -q -E "Job_Name = ${SAMPLE}_links"; then
            progress_log "yellow" "skani for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            new_jobs=$((new_jobs + 1))
            progress_log "magenta" "Submitting links for ${SAMPLE}"
            continue
            sed -i "/${SAMPLE}_links/d" logs/pbs_jobid/links_pbs_jobid.csv
            rm -rf ${OUT_DIR}
            mkdir -p ${OUT_DIR}
        fi
    fi

    qsub <<- EOF | sed -e "s/^/${SAMPLE}_links,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/links_pbs_jobid.csv
#!/bin/env bash

#PBS -N ${SAMPLE}_links
#PBS -j oe
#PBS -o logs/links/${SAMPLE}_links.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

source /home/aubgxs001/miniforge3/bin/activate Links

cd $OUT_DIR

echo $NANOPORE_FQ > ${SAMPLE}_reads.fof
# Run links
LINKS \
    -j 4 \
    -d 500,1000,3000,4000 \
    -f $MEDAKA_FASTA \
    -s ${OUT_DIR}/${SAMPLE}_reads.fof \
    -b ${SAMPLE}.links

if [[ ! -s ${OUT_DIR}/${SAMPLE}.links.scaffolds.fa ]]; then
    echo "Links failed for ${SAMPLE}. Exiting."
    exit 1
else
    echo "Links successful, copying to home"
    cp ${OUT_DIR}/${SAMPLE}.links.scaffolds.fa ${FINAL_DIR}/${SAMPLE}.links.fasta
fi

EOF

progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

done


progress_done "Done submitting links on all samples. Submitted ${new_jobs} new jobs."