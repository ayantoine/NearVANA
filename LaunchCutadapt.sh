#!/bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
source $DATA

echo "------ Get Subsample list ------"
declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	VAR_SAMPLE_FILE="${VARNAME}[3]"
	#echo "${!VAR_SAMPLE_FILE}"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Subsample list ------"

SAMPLE=${SAMPLE_LIST[${STASKID}-1]}

echo "------ Get adaptors ------"
adapter1=$(cut -f1 ${ADAP})
adapter2=$(cut -f2 ${ADAP})
echo "Adapter1: "$adapter1
echo "Adapter2: "$adapter2
echo "------ /Get adaptors ------"

echo "------ Trim adaptors ------"
cutadapt --core=${MULTICPU} -a $adapter1 -A $adapter2 -q 30 -O $((${#adapter1}*85/100)) -m 15 -j 0 -o ${SAMPLE}/${SAMPLE}_${PID}_R1.fastq.split.trim -p ${SAMPLE}/${SAMPLE}_${PID}_R2.fastq.split.trim ${SAMPLE}/${SAMPLE}_${PID}_R1.fastq.split ${SAMPLE}/${SAMPLE}_${PID}_R2.fastq.split
echo "------ /Trim adaptors ------"

touch TrimReads_Ok/${STASKID}.TrimReads.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Cutadapt: "$delta
