#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

nb_jobs=$2

python ${SDIR}/CreateTableFusion.py -j ${nb_jobs} -p ${PID} -m ${META}

touch ${PID}.creationAll.ok
