#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# General Configs
config_walltime="01:00:00"
config_nodes="1"
config_ncpus="4"
config_mem="30gb"
config_queue="large"

# Create output directories
if [[ ! -d logs/skani || ! -d logs/pbs_jobid || -s logs/pbs_jobid/skani_pbs_jobid.csv ]]; then
    mkdir -p logs/skani
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/skani_pbs_jobid.csv
else
    cp logs/pbs_jobid/skani_pbs_jobid.csv logs/pbs_jobid/skani_pbs_jobid.csv.bk
fi

total_checks=$(( ${#FILTERED_SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${FILTERED_SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    # FASTA="${SCRATCHDIR}/fastas/${SAMPLE}.fasta"
    FASTA="${HOMEDIR}/results/assemblies/repeatmasker/links/${SAMPLE}.masked.fasta"
    OUT_DIR="${SCRATCHDIR}/results/skani_links_s40"
    DATABASE="/home/aubgxs001/scratch/databases/skani/database_slow"

    # Check if it has been run before
    if [[ ! -s ${FASTA} ]]; then
        progress_log "red" "skani for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/${SAMPLE}/${SAMPLE}.i.txt" ]]; then
        progress_log "green" "skani for ${SAMPLE} already completed. Skipping."
        continue
    else
        if qstat -f | grep -q -E "Job_Name = skani_${SAMPLE}_links"; then
            progress_log "yellow" "skani for ${SAMPLE} is already running. Skipping."
            continue
        else
            wait_checker
            new_jobs=$((new_jobs + 1))
            progress_log "magenta" "Submitting skani for ${SAMPLE}"
            sed -i "/skani_${SAMPLE}_links,/d" logs/pbs_jobid/skani_pbs_jobid.csv
            mkdir -p "${OUT_DIR}/${SAMPLE}"
        fi
    fi

    qsub > >(sed -e "s/^/skani_${SAMPLE}_links,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/skani_pbs_jobid.csv) 2> >(raise_error) <<- EOF
#!/bin/env bash

#PBS -N skani_${SAMPLE}_links
#PBS -j oe
#PBS -o logs/skani/skani_${SAMPLE}_links.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

source /home/aubgxs001/miniforge3/bin/activate Skani

# Run skani

echo "Running skani for ${SAMPLE}"

cd ${OUT_DIR}

skani search -t 4 \\
    -d ${DATABASE} \\
    -o ${OUT_DIR}/${SAMPLE}/${SAMPLE}.i.txt \\
    -s 40 \\
    ${FASTA}

# skani search -t 4 \\
#     -d ${DATABASE} \\
#     -o ${OUT_DIR}/${SAMPLE}/${SAMPLE}.qi.txt \\
#     --qi ${FASTA}

EOF

progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
done

echo_overwrite_2 "Done submitting skani on all samples. Submitted ${new_jobs} new jobs."