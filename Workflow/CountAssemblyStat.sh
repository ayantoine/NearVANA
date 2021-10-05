#! /bin/bash

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

#N or X or D
TASK=$2

echo "${TASK}"

echo "------ Write stat ------"
echo "python ${SDIR}/CountAssemblyStat.py -0 ${PID}_R0.Substracted.fastq -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -u ${PID}_All.Megahit_unmappedReads.tsv -o ${PID}_Stat_Assembly.tsv"
python ${SDIR}/CountIdentificationStat.py -t ${PID}_Blast${TASK}_results.tab -s ${PID}_Stat_Assembly.tsv -o ${PID}_Stat_identification-${TASK}.tsv
echo "------ /Write stat ------"

echo "------ Create ok tagfile ------"
touch ${PID}.Stat-${TASK}.ok
echo "------ /Create ok tagfile ------"
