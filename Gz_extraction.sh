#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG

zcat ${R1} > ${PID}_R1.fastq
zcat ${R2} > ${PID}_R2.fastq

touch ${PID}.extraction.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Gz_extraction: "$delta > Time01.txt
