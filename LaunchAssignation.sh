#! /bin/bash

set -e

err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

KMER_FILE=$1
WORKDIR=$2
SCRIPTDIR=$3
OKDIR=$4
CONFFILE=$5
ARGFILE=$6

source ${ARGFILE}
source ${CONFFILE}

echo "-1"
echo ${STASKID}

echo "-2"
echo "REALVALUE=$(expr ${STASKID} - 1)"
#Split start at 0 and Task count start at 1. Substract 1 to Task count to match split name
echo "-3"
REALVALUE=$(expr ${STASKID} - 1)
echo $REALVALUE
echo "-4"
DIGITID="000${REALVALUE}"
echo $DIGIT
echo "-5"
DIGITID="${DIGITID: -3}"
echo $DIGIT

echo "-6"
echo "python ${SCRIPTDIR}/MakeAssignation.py -1 ${PID}_R1.fastq -2 ${PID}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID}"
python ${SCRIPTDIR}/MakeAssignation.py -1 ${PID}_R1.fastq -2 ${PID}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID}

echo "-7"
