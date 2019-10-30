#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Build Mapping Reference ------"
bowtie2-build ${SUBS} ${SUBS}
echo "------ /Build Mapping Reference ------"

echo "------ Mapping ------"
bowtie2 --end-to-end --very-sensitive -x ${SUBS} -1 ${PID}_R1.Unsubstracted.fastq -2 ${PID}_R2.Unsubstracted.fastq -U ${PID}_R0.Unsubstracted.fastq -S ${PID}"_"${SUBS}"_"${PID}_R1.Unsubstracted.fastq"_"${PID}_R2.Unsubstracted.fastq"_bwt_ete.sam"
echo "------ /Mapping ------"

echo "------ Mapping extraction ------"
python ${SDIR}/MappingExtraction.py -p ${PID} -i ${PID}"_"${SUBS}"_"${PID}_R1.Unsubstracted.fastq"_"${PID}_R2.Unsubstracted.fastq"_bwt_ete.sam"
echo "------ /Mapping extraction ------"

rm ${PID}"_"${SUBS}"_"${PID}_R1.Unsubstracted.fastq"_"${PID}_R2.Unsubstracted.fastq"_bwt_ete.sam"
rm *.bt2

touch ${PID}.Substraction.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Mapping against PhiX: "$delta > Time06.txt
