#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Correction ------"
spades.py --only-error-correction --disable-gzip-output --pe1-1 ${PID}_R1.Substracted.fastq --pe1-2 ${PID}_R2.Substracted.fastq --pe1-s ${PID}_R0.Substracted.fastq -o ${PID}"_log_Correction"
echo "------ /Correction ------"

cat ${PID}_log_Correction/corrected/${PID}_R1.Substracted.* > ${PID}_R1.Corrected.fastq
cat ${PID}_log_Correction/corrected/${PID}_R2.Substracted.* > ${PID}_R2.Corrected.fastq
cat ${PID}_log_Correction/corrected/${PID}_R0.Substracted.* > ${PID}_R0.Corrected.fastq
cat ${PID}_log_Correction/corrected/${PID}_R_unpaired.* >> ${PID}_R0.Corrected.fastq

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Correction: "$delta > Time07.txt
