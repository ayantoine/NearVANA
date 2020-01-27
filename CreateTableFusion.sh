#! /bin/bash

set -e

err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

nb_jobs=$2

echo "python ${SDIR}/CreateTableFusion.py -j ${nb_jobs} -p ${PID} -m ${META} -l ${VIRMINLEN}"
python ${SDIR}/CreateTableFusion.py -j ${nb_jobs} -p ${PID} -m ${META} -l ${VIRMINLEN}

touch ${PID}.creationAll.ok
