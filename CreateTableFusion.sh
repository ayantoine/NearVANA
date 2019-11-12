#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

task=$2
nb_jobs=$3

python ${SDIR}/CreateTableFusion.py -j ${nb_jobs} -p ${PID} -m ${META}

touch ${PID}.creationAll.ok
