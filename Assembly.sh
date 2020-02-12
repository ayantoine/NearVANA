#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Megahit ------"
megahit --k-list 21,33,55,77,99 -1 VP-1_R1.Corrected.fastq -2 VP-1_R2.Corrected.fastq -r VP-1_R0.Corrected.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${PID}"_log_Assembly-Megahit"
echo "------ /Megahit ------"

mv ${PID}_log_Assembly-Megahit/final.contigs.fa ${PID}_Temp.Megahit_contigs.fa

echo "------ Megahit reverse-mapping ------"
bowtie2-build --threads ${MULTICPU} ${PID}_Temp.Megahit_contigs.fa ${PID}_Temp.Megahit_contigs.fa
bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Corrected.fastq -2 ${PID}_R2.Corrected.fastq -U ${PID}_R0.Corrected.fastq -S reads2contigs.sam
rm *.bt2
python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam
echo "------ /Megahit reverse-mapping ------"

#rm ${PID}_Temp.Megahit_contigs.fa reads2contigs.sam

echo "------ Compress Corrected.fastq ------"
gzip -f ${PID}_R0.Corrected.fastq > ${PID}_R0.Corrected.fastq.gz
gzip -f ${PID}_R1.Corrected.fastq > ${PID}_R1.Corrected.fastq.gz
gzip -f ${PID}_R2.Corrected.fastq > ${PID}_R2.Corrected.fastq.gz
echo "------ /Compress Corrected.fastq ------"

echo "------ FLASH ------"
flash -m 15 -M 300 -O ${PID}_R1.Megahit_unassembled.fastq ${PID}_R2.Megahit_unassembled.fastq -o FLASH
python ${SDIR}/MappingReverseFLASH.py -i FLASH.extendedFrags.fastq -p ${PID}
python ${SDIR}/ConvertFastq2Fasta.py -f FLASH.notCombined_1.fastq -o ${PID}"_R1.FLASH_unassembled.fa"
python ${SDIR}/ConvertFastq2Fasta.py -f FLASH.notCombined_2.fastq -o ${PID}"_R2.FLASH_unassembled.fa"
echo "------ /FLASH ------"

rm FLASH*
rm ${PID}_R1.Megahit_unassembled.fastq ${PID}_R2.Megahit_unassembled.fastq
python ${SDIR}/ConvertFastq2Fasta.py -f ${PID}_R0.Megahit_unassembled.fastq -o ${PID}"_R0.Megahit_unassembled.fa"
rm ${PID}_R0.Megahit_unassembled.fastq

echo "------ Merge Assembly ------"
touch ${PID}_All.fa
cat ${PID}_All.Megahit_contigs.fa >> ${PID}_All.fa
cat ${PID}_All.FLASH_contigs.fa >> ${PID}_All.fa
cat ${PID}_R1.FLASH_unassembled.fa >> ${PID}_All.fa
cat ${PID}_R2.FLASH_unassembled.fa >> ${PID}_All.fa
cat ${PID}_R0.Megahit_unassembled.fa >> ${PID}_All.fa
echo "------ /Merge Assembly ------"

#rm ${PID}_All.Megahit_contigs.fa ${PID}_All.FLASH_contigs.fa ${PID}_R1.FLASH_unassembled.fa
#rm ${PID}_R2.FLASH_unassembled.fa ${PID}_R0.Megahit_unassembled.fa

touch ${PID}.Assembly.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Assembly: "$delta > Time08.txt
