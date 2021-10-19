#! /bin/bash

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

#N or X or D
TASK=$2

echo "${TASK}"

echo "------ Write stat ------"
echo "${SDIR}/CountIdentificationStat.py -t ${PID}_Blast${TASK}_results.tab -s ${PID}_Stat_Assembly.tsv -o ${PID}_Stat_identification-${TASK}.tsv"
python ${SDIR}/CountIdentificationStat.py -t ${PID}_Blast${TASK}_results.tab -s ${PID}_Stat_Assembly.tsv -o ${PID}_Stat_identification-${TASK}.tsv
echo "------ /Write stat ------"

echo "------ Create ok tagfile ------"
touch ${PID}.Stat_Identification-${TASK}.ok
echo "------ /Create ok tagfile ------"
