#!/bin/env bash
set -e

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# General Configs
config_walltime="48:00:00"
config_nodes="1"
config_ncpus="20"
config_mem="300gb"
config_queue="bigmem"

# Create output directories
if [[ ! -d logs/funannotate || ! -d logs/pbs_jobid || -s logs/pbs_jobid/funannotate_pbs_jobid.csv ]]; then
    mkdir -p logs/funannotate
    mkdir -p logs/pbs_jobid
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/funannotate_pbs_jobid.csv
else
    cp logs/pbs_jobid/funannotate_pbs_jobid.csv logs/pbs_jobid/funannotate_pbs_jobid.csv.bk
fi

LOGS_DIR="${HOMEDIR}/scripts/funannotate_links/logs/funannotate"

total_checks=$(( ${#FILTERED_SAMPLES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE in "${FILTERED_SAMPLES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"
    
    GENOME="${HOMEDIR}/results/assemblies/repeatmasker/links/${SAMPLE}.masked.fasta"
    OUT_DIR="${SCRATCHDIR}/results/annotation/funannotate_links/${SAMPLE}"
    RNA_DB="/scratch/aubgxs001/databases/refseq_fungi/rna/fungi.rna.fna"
    PROT_DB="/scratch/aubgxs001/databases/refseq_fungi/protein/fungi.protein.faa"

    mkdir -p ${OUT_DIR}

    # Check if it has been run before
    if [[ ! -s ${GENOME} ]]; then
        progress_log "red" "Funannotate for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/predict/predict_results/fungi_${SAMPLE}.gff3" ]]; then
        progress_log "green" "Funannotate for ${SAMPLE} already completed. Skipping."
        continue
    elif $(tail "${LOGS_DIR}/funannotate_${SAMPLE}_links.out" 2> /dev/null | grep -q "to train Augustus (100 required), exiting"); then
        progress_log "red" "Funannotate for ${SAMPLE} has already been attempted but failed due to insufficient training models. Skipping."
        continue
    elif $(qstat -f | grep -q -E "Job_Name = funannotate_${SAMPLE}_links"); then
        progress_log "yellow" "Funannotate for ${SAMPLE} is already running. Skipping."
        continue
    else
        wait_checker
        progress_log "magenta" "Submitting Funannotate for ${SAMPLE}"
        new_jobs=$((new_jobs + 1))
        sed -i "/funannotate_${SAMPLE}_links,/d" logs/pbs_jobid/funannotate_pbs_jobid.csv
        mkdir -p ${OUT_DIR}
    fi
    
    qsub 2> >(raise_error) > >(sed -e "s/^/funannotate_${SAMPLE}_links,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/funannotate_pbs_jobid.csv) <<- EOF

#!/bin/env bash

#PBS -N funannotate_${SAMPLE}_links
#PBS -j oe
#PBS -o logs/funannotate/funannotate_${SAMPLE}_links.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Funannotate

if [[ ! -s ${OUT_DIR}/${SAMPLE}.cleaned.fasta ]]; then
    funannotate clean --exhaustive --input ${GENOME} --out ${OUT_DIR}/${SAMPLE}.cleaned.fasta
fi

if [[ ! -s ${OUT_DIR}/${SAMPLE}.sorted.fasta ]]; then
    funannotate sort --input ${OUT_DIR}/${SAMPLE}.cleaned.fasta --out ${OUT_DIR}/${SAMPLE}.sorted.fasta --simplify --minlen 500
fi

funannotate predict \
            --cpus 20 \
            --input ${OUT_DIR}/${SAMPLE}.sorted.fasta \
            --out ${OUT_DIR}/predict \
            --species fungi \
            --name ${SAMPLE} \
            --min_training_models 100 \
            --busco_db dikarya \
            --isolate ${SAMPLE} \
            --protein_evidence ${PROT_DB} \
            --transcript_evidence ${RNA_DB} \
            --tmpdir \${workdir}

$(cleanup)
EOF

done


echo_overwrite_2 "Done submitting funannotate on all samples. Submitted ${new_jobs} new jobs."