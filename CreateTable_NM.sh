#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

task=$2
nb_jobs=$3

echo "python ${SDIR}/CreateTable_NM.py -j ${nb_jobs} -p ${PID} -l ${VIRMINLEN}"
python ${SDIR}/CreateTable_NM.py -j ${nb_jobs} -p ${PID} -l ${VIRMINLEN}

touch ${PID}.creation${task}.ok
