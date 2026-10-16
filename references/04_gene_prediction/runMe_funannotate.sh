#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# load_env Flye

# General Configs
config_walltime="24:00:00"
config_nodes="1"
config_ncpus="18"
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

# Read the list of samples from the CSV file, skipping the header and make spaces in the name underscores
readarray -t REFERENCES < <(awk -F',' 'NR > 1 {print $1}' /home/aubgxs001/scratch/comparative_fungi/references/reference_genomes.csv)
REFERENCES=("${REFERENCES[@]// /_}")

total_checks=$(( ${#REFERENCES[@]}))
echo "Total checks to consider: ${total_checks}"

for SAMPLE_FULL in "${REFERENCES[@]}"; do
    sample_number=$(( sample_number + 1 ))
    current_check=$(( current_check + 1 ))
    progress_bar "$current_check" "$total_checks" "$start_time" "$sample_number" "New jobs: $new_jobs"

    GENOME="${SCRATCHDIR}/references/${SAMPLE_FULL}/${SAMPLE_FULL}_genomic.fna"
    SAMPLE=$(echo "${SAMPLE_FULL}" | tr '_' '-' | awk -F'-' '{print substr($1,1,1) $2}')
    OUT_DIR="${SCRATCHDIR}/results/annotation/funannotate/${SAMPLE}"
    RNA_DB="/scratch/aubgxs001/databases/refseq_fungi/rna/fungi.rna.fna"
    PROT_DB="/scratch/aubgxs001/databases/refseq_fungi/protein/fungi.protein.faa"
    # Abreviate by takin the first letter of the first word and append to the second word, with no separator.
    # SAMPLE_ABV=$(echo "${SAMPLE}" | tr '_' '-' | awk -F'-' '{print substr($1,1,1) "-" substr($2,1,3)}')

    mkdir -p ${OUT_DIR}
    # Check if it has been run before
    if [[ ! -s ${GENOME} ]]; then
        progress_log "red" "Funannotate for ${SAMPLE} cannot run because one or more of the input files does not exist. Skipping."
        continue
    elif [[ -s "${OUT_DIR}/predict/predict_results/fungi_${SAMPLE}.gff3" ]]; then
        progress_log "green" "Funannotate for ${SAMPLE} already completed. Skipping."
        continue
    # elif $(tail "${HOMEDIR}/scripts/logs/funannotate/funannotate_${SAMPLE}.out" | grep -q "to train Augustus (100 required), exiting"); then
    #     progress_log "red" "Funannotate for ${SAMPLE} has already been attempted but failed due to insufficient training models. Skipping."
    #     continue
    elif qstat -f | grep -q -E "Job_Name = funannotate_${SAMPLE}"; then
        progress_log "yellow" "Funannotate for ${SAMPLE} is already running. Skipping."
        continue
    else
        wait_checker
        progress_log "magenta" "Submitting Funannotate for ${SAMPLE}"
        new_jobs=$((new_jobs + 1))
        sed -i "/funannotate_${SAMPLE},/d" logs/pbs_jobid/funannotate_pbs_jobid.csv
        # rm -rf ${OUT_DIR}
        # mkdir -p ${OUT_DIR}
    fi

    qsub <<- EOF | sed -e "s/^/funannotate_${SAMPLE},/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/funannotate_pbs_jobid.csv
#!/bin/env bash

#PBS -N funannotate_${SAMPLE}
#PBS -j oe
#PBS -o logs/funannotate/funannotate_${SAMPLE}.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Funannotate

cut -d ' ' -f1 ${GENOME} > ${OUT_DIR}/${SAMPLE}.space.fna

seqkit seq ${OUT_DIR}/${SAMPLE}.space.fna -w 60 > ${OUT_DIR}/${SAMPLE}.contigs.fna

funannotate clean --exhaustive --input ${OUT_DIR}/${SAMPLE}.contigs.fna --out ${OUT_DIR}/${SAMPLE}.cleaned.fasta

# funannotate sort --input ${OUT_DIR}/${SAMPLE}.cleaned.fasta --out ${OUT_DIR}/${SAMPLE}.sorted.fasta --simplify --minlen 500

mkdir ${OUT_DIR}/predict

funannotate predict \
            --cpus 20 \
            --input ${OUT_DIR}/${SAMPLE}.cleaned.fasta \
            --out ${OUT_DIR}/predict \
            --species fungi \
            --name ${SAMPLE//_/-} \
            --min_training_models 100 \
            --busco_db dikarya \
            --isolate ${SAMPLE} \
            --header_length 400 \
            --protein_evidence ${PROT_DB} \
            --transcript_evidence ${RNA_DB} \
            --tmpdir \${workdir}

$(cleanup)

EOF
done


echo_overwrite_2 "Done submitting funannotate on all samples. Submitted ${new_jobs} new jobs."