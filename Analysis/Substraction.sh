#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Analysis

echo "------ Build Mapping Reference ------"
bowtie2-build --threads ${MULTICPU} ${SUBS} ${SUBS}
echo "------ /Build Mapping Reference ------"

echo "------ Mapping ------"
bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${SUBS} -1 ${PID}_R1.Unsubstracted.fastq -2 ${PID}_R2.Unsubstracted.fastq -U ${PID}_R0.Unsubstracted.fastq -S ${PID}"_"${SUBS}"_"${PID}_R1.Unsubstracted.fastq"_"${PID}_R2.Unsubstracted.fastq"_bwt_ete.sam"
echo "------ /Mapping ------"

echo "------ Mapping extraction ------"
python ${SDIR}/MappingExclusion.py -p ${PID} -i ${PID}"_"${SUBS}"_"${PID}_R1.Unsubstracted.fastq"_"${PID}_R2.Unsubstracted.fastq"_bwt_ete.sam"
echo "------ /Mapping extraction ------"

#rm ${PID}"_"${SUBS}"_"${PID}_R1.Unsubstracted.fastq"_"${PID}_R2.Unsubstracted.fastq"_bwt_ete.sam"
#rm *.bt2

echo "------ Store data ------"
gzip -f ${PID}_R0.PhiX.fastq > ${PID}_R0.PhiX.fastq.gz
gzip -f ${PID}_R1.PhiX.fastq > ${PID}_R1.PhiX.fastq.gz
gzip -f ${PID}_R2.PhiX.fastq > ${PID}_R2.PhiX.fastq.gz
echo "------ Store data ------"

touch ${PID}.Substraction.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Mapping against PhiX: "$delta > Time06.txt
