#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"

echo "------ Megahit ------"
if [ ! -f ${PID}_Temp.Megahit_contigs.fa ]; then
    time megahit --k-list 21,33,55,77,99 -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -r ${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${PID}"_log_Assembly-Megahit"
    mv ${PID}_log_Assembly-Megahit/final.contigs.fa ${PID}_Temp.Megahit_contigs.fa
fi
echo "------ /Megahit ------"

echo "------ Megahit reverse-mapping ------"
if [ ! -f ${PID}_All.Megahit_contigs.fa ]; then
    bowtie2-build --threads ${MULTICPU} ${PID}_Temp.Megahit_contigs.fa ${PID}_Temp.Megahit_contigs.fa
    if [ -s "${PID}_R0.Substracted.fastq" ]; then
	    bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -U ${PID}_R0.Substracted.fastq -S reads2contigs.sam
    else
	    bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -S reads2contigs.sam
    fi
    rm *.bt2
    echo "python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}"
    python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}
    rm ${PID}_Temp.Megahit_contigs.fa reads2contigs.sam
fi
echo "------ /Megahit reverse-mapping ------"

if [ "$USE_MULTIPLEX" = true ] ; then
	if [ ! -f ${PID}.Stat_Assembly.ok ]; then
		echo "------ Write stat ------"
		echo "$SCALL $SPARAM_HEAVY $SRENAME ${PID}_${TASK}Stat -e Stat_Assembly.e -o Stat_Assembly.o ${SDIR}/CountAssemblyStat.sh $ARG"
		$SCALL $SPARAM_HEAVY $SRENAME ${PID}_${TASK}Stat -e Stat_Assembly.e -o Stat_Assembly.o ${SDIR}/CountAssemblyStat.sh $ARG
		while [ ! -e ${PID}.Stat_Assembly.ok ]; do sleep 60 ; done
		echo "------ /Write stat ------"
	fi
fi

echo "------ Compress Corrected.fastq ------"
gzip -f ${PID}_R0.Substracted.fastq
gzip -f ${PID}_R1.Substracted.fastq
gzip -f ${PID}_R2.Substracted.fastq
echo "------ /Compress Corrected.fastq ------"

echo "------ Merge Assembly ------"
mv ${PID}_All.Megahit_contigs.fa ${PID}_All.fa
echo "------ /Merge Assembly ------"

echo "------ Compress Megahit output ------"
gzip -f ${PID}_All.Megahit_reverseAssembly.tsv
gzip -f ${PID}_All.Megahit_ambigousReads.tsv
gzip -f ${PID}_All.Megahit_rejectedContigs.fa
gzip -f ${PID}_All.Megahit_unmappedReads.tsv
echo "------ /Compress Megahit output ------"

touch ${PID}.Assembly.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Assembly: "$delta > Time07.txt
