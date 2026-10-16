#!/bin/env bash
set -e

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=48:00:00
config_ncpus=35
config_mem=80gb
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
    GENOME="${SCRATCHDIR}/references/${SAMPLE_FULL}/${SAMPLE_FULL}_genomic.fna"
    OUT_DIR="${SCRATCHDIR}/results/annotation/funannotate/${SAMPLE}"

    mkdir -p ${OUT_DIR}

    # Check if it has been run before
    if [[ ! -s ${GENOME} ]]; then
        progress_log "red" "Funannotate for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ ! -s "${OUT_DIR}/predict/predict_results/fungi_${SAMPLE}.gff3" ]]; then
        progress_log "red" "Funannotate Interpro for ${SAMPLE} cannot be run because annotation is not available. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/predict/annotate_misc/iprscan.xml" ]]; then
        progress_log "green" "Funannotate Interpro for ${SAMPLE} has already completed. Skipping."
        continue
    # elif qstat -f | grep -q -E "Job_Name = funannotate_iprscan_${SAMPLE}"; then
    #     progress_log "yellow" "Funannotate Interpro for ${SAMPLE} is already running. Skipping."
    #     continue
    else
        wait_checker
        progress_log "magenta" "Submitting Funannotate for ${SAMPLE}"
        new_jobs=$((new_jobs + 1))
        sed -i "/funannotate_iprscan_${SAMPLE},/d" logs/pbs_jobid/funannotate_pbs_jobid.csv
    fi

    qsub <<- EOF | sed -e "s/^/funannotate_iprscan_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/funannotate_pbs_jobid.csv
#!/bin/env bash

#PBS -N funannotate_iprscan_${SAMPLE}
#PBS -j oe
#PBS -o logs/funannotate/funannotate_iprscan_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Funannotate

funannotate iprscan -i ${OUT_DIR}/predict \
            --cpus 10 \
            --m local \
            --iprscan_path "/home/aubgxs001/scratch/tools/interproscan-5.77-108.0/interproscan.sh"

$(cleanup)
EOF

done
echo_overwrite_2 "Done submitting funannotate on all samples. Submitted ${new_jobs} new jobs."