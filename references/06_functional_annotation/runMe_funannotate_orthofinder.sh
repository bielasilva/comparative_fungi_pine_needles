#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

# load_env Flye

# General Configs
config_walltime="12:00:00"
config_nodes="1"
config_ncpus="32"
config_mem="120gb"
config_queue="large"

# Create output directories
if [[ ! -d logs/orthofinder || ! -d logs/pbs_jobid || -s logs/pbs_jobid/orthofinder_pbs_jobid.csv ]]; then
    mkdir -p logs/orthofinder
    mkdir -p logs/pbs_jobid
    mkdir -p submission/orthofinder
    echo -e "JobName,JobID,SubmissionTime" > logs/pbs_jobid/orthofinder_pbs_jobid.csv
else
    cp logs/pbs_jobid/orthofinder_pbs_jobid.csv logs/pbs_jobid/orthofinder_pbs_jobid.csv.bk
fi


# mkdir -p ${SCRATCHDIR}/results/orthofinder/proteins
# ln -s ${SCRATCHDIR}/results/annotation/funannotate/*/predict/annotate_results/*.proteins.fa ${SCRATCHDIR}/results/orthofinder/proteins/

rm -rf /scratch/aubgxs001/comparative_fungi/results/orthofinder/orthology/

# # Check if it has been run before
# if [[ ! -s "${PROT}" ]]; then
#     progress_log "red" "Proteins for ${SAMPLE} are not available. Skipping."
#     # rm -rf ${FUN_DIR}
#     continue
# elif qstat -f | grep -q -E "Job_Name = orthofinder_${SAMPLE}"; then
#     progress_log "yellow" "orthofinder Annotate for ${SAMPLE} is already running. Skipping."
#     continue
# elif qstat -f | grep -q -E "Job_Name = orthofinder_${SAMPLE}"; then
#     progress_log "yellow" "orthofinder for ${SAMPLE} is still running. Skipping."
#     continue
# else
#     wait_checker
#     new_jobs=$((new_jobs + 1))
#     progress_log "magenta" "Submitting orthofinder for ${SAMPLE}"
#     continue
#     sed -i "/orthofinder_${SAMPLE},/d" logs/pbs_jobid/orthofinder_pbs_jobid.csv
# fi

qsub <<- EOF | sed -e "s/^/orthofinder,/" -e "s/$/,$(date +'%Y-%m-%d %H:%M:%S')/" >> logs/pbs_jobid/orthofinder_pbs_jobid.csv
#!/bin/env bash

#PBS -N orthofinder
#PBS -j oe
#PBS -o logs/orthofinder/orthofinder.out
#PBS -l walltime=${config_walltime}
#PBS -l select=ncpus=${config_ncpus}:mpiprocs=${config_ncpus}:mem=${config_mem}
#PBS -q ${config_queue}

$(log_asax)

source /home/aubgxs001/miniforge3/bin/activate Orthofinder

orthofinder \
    -f ${SCRATCHDIR}/results/orthofinder/proteins \
    -t 32 \
    -a 32 \
    -o ${SCRATCHDIR}/results/orthofinder/orthology

EOF



echo_overwrite_2 "Done submitting orthofinder on all samples. Submitted ${new_jobs} new jobs."