#! /bin/bash

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

#N or X or D
TASK=$2

echo "${TASK}"

echo "------ Xlsx conversion ------"
echo "perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_Blast${TASK}_results.tab ${PID}_Blast${TASK}_results.xlsx $((${#PID}+7))"
perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_Blast${TASK}_results.tab ${PID}_Blast${TASK}_results.xlsx $((${#PID}+7))
echo "------ /Xlsx conversion ------"


echo "------ Create ok tagfile ------"
touch ${PID}.Tab2Xls${TASK}.ok
echo "------ /Create ok tagfile ------"
