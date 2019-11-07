#! /bin/bash

KMER_FILE=$1
WORKDIR=$2
SCRIPTDIR=$3
OKDIR=$4
CONFFILE=$5
ARGFILE=$6

source ${ARGFILE}
source ${CONFFILE}

REALVALUE=$(expr ${STASKID} - 1) #Split start at 0 and Task count start at 1. Substract 1 to Task count to match split name
DIGITID="000${REALVALUE}"
DIGITID="${DIGITID: -3}"


python ${SCRIPTDIR}/MakeAssignation.py -1 ${PID}_R1.fastq -2 ${PID}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID}