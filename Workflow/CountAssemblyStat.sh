#! /bin/bash

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

echo "------ Write stat ------"
echo "python ${SDIR}/CountAssemblyStat.py -0 ${PID}_R0.Substracted.fastq -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -u ${PID}_All.Megahit_unmappedReads.tsv -o ${PID}_Stat_Assembly.tsv -d ${PID}_Demultiplexing_Hyper_Distribution.tsv"
python ${SDIR}/CountAssemblyStat.py -0 ${PID}_R0.Substracted.fastq -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -u ${PID}_All.Megahit_unmappedReads.tsv -o ${PID}_Stat_Assembly.tsv -d ${PID}_Demultiplexing_Hyper_Distribution.tsv
echo "------ /Write stat ------"

echo "------ Create ok tagfile ------"
touch ${PID}.Stat_Assembly.ok
echo "------ /Create ok tagfile ------"
