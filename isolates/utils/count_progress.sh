#!/bin/env bash

source /home/aubgxs001/comparative_fungi/scripts/sourceMe_configs.sh

flye_count=0
medaka01_count=0
medaka02_count=0
links_count=0
fun_predict_count=0
antismash_count=0
signalp_count=0
dtm2_count=0
effectorp_count=0
interpro_count=0
fun_annotate_count=0

for SAMPLE in "${FILTERED_SAMPLES[@]}"; do

    FUN_DIR="${SCRATCHDIR}/results/annotation/funannotate_links/${SAMPLE}"

    if [[ -s "${HOMEDIR}/results/assemblies/flye/${SAMPLE}.flye.fasta" ]]; then
        ((flye_count += 1))
    fi

    if [[ -s "${HOMEDIR}/results/assemblies/medaka/01/${SAMPLE}.medaka01.fasta" ]]; then
        ((medaka01_count += 1))
    fi

    if [[ -s "${HOMEDIR}/results/assemblies/medaka/02/${SAMPLE}.medaka02.fasta" ]]; then
        ((medaka02_count += 1))
    fi

    if [[ -s "${HOMEDIR}/results/assemblies/links/${SAMPLE}.links.fasta" ]]; then
        ((links_count += 1))
    fi
    
    if [[ -s "${FUN_DIR}/predict/predict_results/fungi_${SAMPLE}.gff3" ]]; then
        ((fun_predict_count += 1))
    fi

    if [[ -s "${FUN_DIR}/antismash/fungi_${SAMPLE}.gbk" ]]; then
        ((antismash_count += 1))
    fi
    
    if [[ -s "${FUN_DIR}/signalp/prediction_results.txt" ]]; then
        ((signalp_count += 1))
    fi

    if [[ -s "${FUN_DIR}/deeptmhmm/filtered_signalp.fasta" ]]; then
        ((dtm2_count += 1))
    fi

    if [[ -s "${FUN_DIR}/effectorp/${SAMPLE}.tsv" ]]; then
        ((effectorp_count += 1))
    fi

    if [[ -s "${FUN_DIR}/predict/annotate_misc/iprscan.xml" ]]; then
        ((interpro_count += 1))
    fi

    if [[ -s "${FUN_DIR}/predict/annotate_results/${SAMPLE}_${SAMPLE}.gff3" ]]; then
        ((fun_annotate_count += 1))
    fi

done

echo "The count for flye is ${flye_count}"
echo "The count for medaka01 is ${medaka01_count}"
echo "The count for medaka02 is ${medaka02_count}"
echo "The count for links is ${links_count}"
echo "The count for fun_predict is ${fun_predict_count}"
echo "The count for antismash is ${antismash_count}"
echo "The count for signalp is ${signalp_count}"
echo "The count for dtm2 is ${dtm2_count}"
echo "The count for effectorp is ${effectorp_count}"
echo "The count for interpro is ${interpro_count}"
echo "The count for fun_annotate is ${fun_annotate_count}"