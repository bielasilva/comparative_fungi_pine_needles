#!/bin/env bash
set -e

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# Define general configs
config_walltime=48:00:00
config_ncpus=20
config_mem=50gb
config_queue=large

# Create output directories
if [[ ! -d logs/funannotate || ! -d logs/pbs_jobid || -s logs/pbs_jobid/funannotate_pbs_jobid.csv ]]; then
    mkdir -p logs/funannotate
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/funannotate_pbs_jobid.csv
else
    cp logs/pbs_jobid/funannotate_pbs_jobid.csv logs/pbs_jobid/funannotate_pbs_jobid.csv.bk
fi

total_checks=$(( ${#FILTERED_SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${FILTERED_SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    GENOME="${HOMEDIR}/results/assemblies/repeatmasker/links/${SAMPLE}.masked.fasta"
    FUN_DIR="${SCRATCHDIR}/results/annotation/funannotate_links/${SAMPLE}"

    # Check if it has been run before
    if [[ ! -s ${GENOME} ]]; then
        progress_log "red" "Funannotate Interpro for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ ! -s "${FUN_DIR}/predict/predict_results/fungi_${SAMPLE}.gff3" ]]; then
        progress_log "red" "Funannotate Interpro for ${SAMPLE} cannot be run because annotation is not available. Skipping."
        continue
    elif [[ -s "${FUN_DIR}/predict/annotate_misc/iprscan.xml" ]]; then
        progress_log "green" "Funannotate Interpro for ${SAMPLE} has already completed. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = iprscan_${SAMPLE}_links"; then
        progress_log "yellow" "Funannotate Interpro for ${SAMPLE} is already running. Skipping."
        continue
    else
        wait_checker
        new_jobs=$((new_jobs + 1))
        progress_log "magenta" "Submitting Funannotate Interpro for ${SAMPLE}"
        sed -i "/iprscan_${SAMPLE}_links,/d" logs/pbs_jobid/funannotate_pbs_jobid.csv
    fi

    qsub  2> >(raise_error) > >(sed -e "s/^/iprscan_${SAMPLE}_links,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/funannotate_pbs_jobid.csv) <<- EOF
#!/bin/env bash

#PBS -N iprscan_${SAMPLE}_links
#PBS -j oe
#PBS -o logs/funannotate/iprscan_${SAMPLE}_links.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Funannotate

funannotate iprscan -i ${FUN_DIR}/predict \
            --cpus 8 \
            --m local \
            --iprscan_path "/home/aubgxs001/scratch/tools/interproscan-5.77-108.0/interproscan.sh"

$(cleanup)
EOF

done
echo_overwrite_2 "Done submitting funannotate on all samples. Submitted ${new_jobs} new jobs."