#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ SPAdes ------"
spades.py -k 11,21,33,55,77,99 --only-assembler --pe1-1 ${PID}"_R1.Corrected.fastq" --pe1-2 ${PID}"_R2.Corrected.fastq" --pe1-s ${PID}"_R0.Corrected.fastq" -o ${PID}"_log_Assembly-SPAdes"
echo "------ /SPAdes ------"

mv ${PID}_log_Assembly-SPAdes/contigs.fasta ${PID}_Temp.SPAdes_contigs.fa
rm -r ${PID}_log_Assembly-SPAdes/K* ${PID}_log_Assembly-SPAdes/misc ${PID}_log_Assembly-SPAdes/tmp
rm ${PID}_log_Assembly-SPAdes/*fast* ${PID}_log_Assembly-SPAdes/*paths

echo "------ SPAdes reverse-mapping ------"
bowtie2-build ${PID}_Temp.SPAdes_contigs.fa ${PID}_Temp.SPAdes_contigs.fa
bowtie2 --end-to-end --very-sensitive -x ${PID}_Temp.SPAdes_contigs.fa -1 ${PID}_R1.Corrected.fastq -2 ${PID}_R2.Corrected.fastq -U ${PID}_R0.Corrected.fastq -S reads2contigs.sam
rm *.bt2
python ${SDIR}/MappingReverseSPAdes.py -p ${PID} -i reads2contigs.sam
echo "------ /SPAdes reverse-mapping ------"

rm ${PID}_Temp.SPAdes_contigs.fa reads2contigs.sam

echo "------ Compress Corrected.fastq ------"
gzip ${PID}_R0.Corrected.fastq > ${PID}_R0.Corrected.fastq.gz
gzip ${PID}_R1.Corrected.fastq > ${PID}_R1.Corrected.fastq.gz
gzip ${PID}_R2.Corrected.fastq > ${PID}_R2.Corrected.fastq.gz
echo "------ /Compress Corrected.fastq ------"

echo "------ FLASH ------"
flash -m 15 -M 300 -O ${PID}_R1.SPAdes_unassembled.fastq ${PID}_R2.SPAdes_unassembled.fastq -o FLASH
python ${SDIR}/MappingReverseFLASH.py -i FLASH.extendedFrags.fastq -p ${PID}
python ${SDIR}/ConvertFastq2Fasta.py -f FLASH.notCombined_1.fastq -o ${PID}"_R1.FLASH_unassembled.fa"
python ${SDIR}/ConvertFastq2Fasta.py -f FLASH.notCombined_2.fastq -o ${PID}"_R2.FLASH_unassembled.fa"
echo "------ /FLASH ------"

rm FLASH*
rm ${PID}_R1.SPAdes_unassembled.fastq ${PID}_R2.SPAdes_unassembled.fastq
python ${SDIR}/ConvertFastq2Fasta.py -f ${PID}_R0.SPAdes_unassembled.fastq -o ${PID}"_R0.SPAdes_unassembled.fa"
rm ${PID}_R0.SPAdes_unassembled.fastq

echo "------ Merge Assembly ------"
touch ${PID}_All.fa
cat ${PID}_All.SPAdes_contigs.fa >> ${PID}_All.fa
cat ${PID}_All.FLASH_contigs.fa >> ${PID}_All.fa
cat ${PID}_R1.FLASH_unassembled.fa >> ${PID}_All.fa
cat ${PID}_R2.FLASH_unassembled.fa >> ${PID}_All.fa
cat ${PID}_R0.SPAdes_unassembled.fa >> ${PID}_All.fa
echo "------ /Merge Assembly ------"

rm ${PID}_All.SPAdes_contigs.fa ${PID}_All.FLASH_contigs.fa ${PID}_R1.FLASH_unassembled.fa
rm ${PID}_R2.FLASH_unassembled.fa ${PID}_R0.SPAdes_unassembled.fastq

touch ${PID}.Assembly.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Assembly: "$delta > Time08.txt
