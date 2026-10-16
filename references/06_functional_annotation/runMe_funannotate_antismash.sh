#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=24:00:00
config_ncpus=5
config_mem=120gb
config_queue=large

# Create output directories
if [[ ! -d logs/funannotate || ! -d logs/pbs_jobid || -s logs/pbs_jobid/funannotate_pbs_jobid.csv ]]; then
    mkdir -p logs/funannotate
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/funannotate_pbs_jobid.csv
else
    cp logs/pbs_jobid/funannotate_pbs_jobid.csv logs/pbs_jobid/funannotate_pbs_jobid.csv.bk
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
    OUT_DIR="${FUN_DIR}/antismash/"
    GBK_FILE="${FUN_DIR}/predict/predict_results/fungi_${SAMPLE}.gbk"
    GFF_FILE="${FUN_DIR}/predict/predict_results/fungi_${SAMPLE}.gff3"

    mkdir -p ${OUT_DIR}

    # Check if it has been run before
    if [[ ! -s "${GFF_FILE}" ]]; then
        progress_log "red" "Funannotate Antismash for ${SAMPLE} cannot be run because annotation is not available. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/fungi_${SAMPLE}.gbk" ]]; then
        progress_log "green" "Funannotate Antismash for ${SAMPLE} is complete. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = funannotate_antismash_${SAMPLE}"; then
        progress_log "yellow" "Funannotate Antismash for ${SAMPLE} is already running. Skipping."
        continue
    else
        wait_checker
        progress_log "magenta" "Submitting Funannotate for ${SAMPLE}"
        continue
        new_jobs=$((new_jobs + 1))
        sed -i "/funannotate_antismash_${SAMPLE},/d" logs/pbs_jobid/funannotate_pbs_jobid.csv
    fi

    qsub <<- EOF | sed -e "s/^/funannotate_antismash_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/funannotate_pbs_jobid.csv
#!/bin/env bash

#PBS -N funannotate_antismash_${SAMPLE}
#PBS -j oe
#PBS -o logs/funannotate/funannotate_antismash_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Antismash
export HMMER_NCPU=20

antismash ${GBK_FILE} \
    --cpus 20 \
    --taxon fungi \
    --genefinding-tool none \
    --no-zip-output \
    --output-dir \${workdir}

if [ $? -eq 0 ]; then
    echo "Funannotate Antismash for ${SAMPLE} completed successfully."
    cp \${workdir}/fungi_${SAMPLE}.{json,gbk} ${OUT_DIR}/
    if [[ -s ${OUT_DIR}/fungi_${SAMPLE}.json ]]; then
        $(cleanup)
    else
        echo "Copying wasn't successful"
    fi
else
    echo "Funannotate Antismash for ${SAMPLE} Failed with exit code $?"
fi

$(cleanup)
EOF

done


echo_overwrite_2 "Done submitting funannotate on all samples. Submitted ${new_jobs} new jobs."