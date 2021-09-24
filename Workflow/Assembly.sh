#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

echo "------ Megahit ------"
time megahit --k-list 21,33,55,77,99 -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -r ${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${PID}"_log_Assembly-Megahit"
echo "------ /Megahit ------"

mv ${PID}_log_Assembly-Megahit/final.contigs.fa ${PID}_Temp.Megahit_contigs.fa

echo "------ Megahit reverse-mapping ------"
bowtie2-build --threads ${MULTICPU} ${PID}_Temp.Megahit_contigs.fa ${PID}_Temp.Megahit_contigs.fa
if [ -s "${PID}_R0.Substracted.fastq" ]; then
	bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -U ${PID}_R0.Substracted.fastq -S reads2contigs.sam
else
	bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -S reads2contigs.sam
fi
rm *.bt2
echo "python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}"
python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}
echo "------ /Megahit reverse-mapping ------"

rm ${PID}_Temp.Megahit_contigs.fa reads2contigs.sam

echo "------ Compress Corrected.fastq ------"
gzip -f ${PID}_R0.Substracted.fastq > ${PID}_R0.Substracted.fastq.gz
gzip -f ${PID}_R1.Substracted.fastq > ${PID}_R1.Substracted.fastq.gz
gzip -f ${PID}_R2.Substracted.fastq > ${PID}_R2.Substracted.fastq.gz
echo "------ /Compress Corrected.fastq ------"

echo "------ Merge Assembly ------"
mv ${PID}_All.Megahit_contigs.fa ${PID}_All.fa
echo "------ /Merge Assembly ------"

echo "------ Write stat ------"
echo "python ${SDIR}/CountAssemblyStat.py -0 ${PID}_R0.Substracted.fastq -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -u ${PID}_All.Megahit_unmappedReads.tsv -o ${PID}_Stat_Assembly.tsv"
python ${SDIR}/CountAssemblyStat.py -0 ${PID}_R0.Substracted.fastq -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -u ${PID}_All.Megahit_unmappedReads.tsv -o ${PID}_Stat_Assembly.tsv
echo "------ /Write stat ------"

echo "------ Compress Megahit output ------"
gzip -f ${PID}_All.Megahit_reverseAssembly.tsv > ${PID}_All.Megahit_reverseAssembly.tsv.gz
gzip -f ${PID}_All.Megahit_ambigousReads.tsv > ${PID}_All.Megahit_ambigousReads.tsv.gz
gzip -f ${PID}_All.Megahit_rejectedContigs.fa > ${PID}_All.Megahit_rejectedContigs.fa.gz
gzip -f ${PID}_All.Megahit_unmappedReads.tsv > ${PID}_All.Megahit_unmappedReads.tsv.gz
echo "------ /Compress Megahit output ------"

touch ${PID}.Assembly.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Assembly: "$delta > Time07.txt
