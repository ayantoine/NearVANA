#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $DATA

for VARNAME in "${PLATE[@]}"; do
	VAR_R1_FILE="${VARNAME}[0]"
	VAR_R2_FILE="${VARNAME}[1]"
	
	zcat ${!VAR_R1_FILE} > ${PID}_${VARNAME}_R1.fastq
	zcat ${!VAR_R2_FILE} > ${PID}_${VARNAME}_R2.fastq
	
done

touch ${PID}.extraction.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Gz_extraction: "$delta > Time01.txt
