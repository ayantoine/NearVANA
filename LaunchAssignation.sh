#! /bin/bash

SEQ_NUMBER=$1
KMER_FILE=$2
WORKDIR=$3
SCRIPTDIR=$4
OKDIR=$5
CONFFILE=$6
ARGFILE=$7

source ${ARGFILE}
source ${CONFFILE}

REALVALUE=$(expr ${STASKID} - 1) #Split start at 0 and Task count start at 1. Substract 1 to Task count to match split name
DIGITID="0000${REALVALUE}"
DIGITID="${DIGITID: -4}"

echo "python ${SCRIPTDIR}/MakeAssignation.py -q ${SEQ_NUMBER} -1 ${PID}_R1.fastq -2 ${PID}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID}"
python ${SCRIPTDIR}/MakeAssignation.py -q ${SEQ_NUMBER} -1 ${PID}_R1.fastq -2 ${PID}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}
