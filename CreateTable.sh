#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

nb_jobs=$2

echo "python ${SDIR}/CreateTable.py -j ${nb_jobs} -p ${PID} -d ${DATA} -l ${VIRMINLEN}"
python ${SDIR}/CreateTable.py -j ${nb_jobs} -p ${PID} -d ${DATA} -l ${VIRMINLEN}

touch ${PID}.creation${task}.ok
