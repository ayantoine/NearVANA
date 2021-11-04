#! /bin/bash

ARG=$1
TASK=$2
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow


echo "------ Write stat ------"
echo "python ${SDIR}/CountIdentificationStat.py -t ${PID}_Blast${TASK}_results.tab -s ${PID}_Stat_Assembly.tsv -o ${PID}_Stat_Identificaion-${TASK}.tsv"
python ${SDIR}/CountIdentificationStat.py -t ${PID}_Blast${TASK}_results.tab -s ${PID}_Stat_Assembly.tsv -o ${PID}_Stat_Identificaion-${TASK}.tsv
echo "------ /Write stat ------"

echo "------ Create ok tagfile ------"
touch ${PID}.Stat_Identification-${TASK}.ok
echo "------ /Create ok tagfile ------"
