#!/bin/env bash
set -e

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# General Configs
config_walltime="12:00:00"
config_nodes="1"
config_ncpus="20"
config_mem="60gb"
config_queue="large"

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
    ANTISMASH="${FUN_DIR}/antismash/fungi_${SAMPLE}.gbk"
    INTERPRO="${FUN_DIR}/predict/annotate_misc/iprscan.xml"
    SIGNALP="${FUN_DIR}/signalp/prediction_results.txt"
    RNA_DB="/scratch/aubgxs001/databases/refseq_fungi/rna/fungi.rna.fna"
    PROT_DB="/scratch/aubgxs001/databases/refseq_fungi/protein/fungi.protein.faa"

    # Check if it has been run before
    if [[ ! -s ${GENOME} ]]; then
        progress_log "red" "Funannotate for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ ! -s $ANTISMASH ]]; then
        progress_log "red" "Antismash results for ${SAMPLE} are not available. Skipping."
        continue
    elif [[ ! -s $INTERPRO ]]; then
        progress_log "red" "Interproscan results for ${SAMPLE} are not available. Skipping."
        continue
    elif [[ ! -s $SIGNALP ]]; then
        progress_log "red" "SignalP results for ${SAMPLE} are not available. Skipping."
        continue
    elif [[ -s "${FUN_DIR}/predict/annotate_results/${SAMPLE}_${SAMPLE}.gff3" ]]; then
        progress_log "green" "Funannotate results for ${SAMPLE} already completed. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = funannotate_annotate_${SAMPLE}"; then
        progress_log "yellow" "Funannotate Annotate for ${SAMPLE} is already running. Skipping."
        continue
    elif qstat -f | grep -q -E "Job_Name = funannotate_${SAMPLE}"; then
        progress_log "yellow" "Funannotate for ${SAMPLE} is still running. Skipping."
        continue
    else
        wait_checker
        progress_log "magenta" "Submitting Funannotate for ${SAMPLE}"
        new_jobs=$((new_jobs + 1))
        sed -i "/funannotate_annotate_${SAMPLE},/d" logs/pbs_jobid/funannotate_pbs_jobid.csv
    fi

    qsub <<- EOF | sed -e "s/^/funannotate_annotate_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/funannotate_pbs_jobid.csv
#!/bin/env bash

#PBS -N funannotate_annotate_${SAMPLE}
#PBS -j oe
#PBS -o logs/funannotate/funannotate_annotate_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Funannotate

funannotate annotate \
            --cpus 20 \
            --input ${FUN_DIR}/predict \
            --species ${SAMPLE} \
            --busco_db dikarya \
            --antismash ${ANTISMASH} \
            --iprscan ${INTERPRO} \
            --signalp ${SIGNALP} \
            --header_length 40

$(cleanup)
EOF
done


echo_overwrite_2 "Done submitting funannotate on all samples. Submitted ${new_jobs} new jobs."