#! /bin/bash

SEQ_NUMBER=$1
KMER_FILE=$2
WORKDIR=$3
SCRIPTDIR=$4
OKDIR=$5
CONFFILE=$6
ARGFILE=$7
VARNAME=$8

source ${ARGFILE}
source ${CONFFILE}

REALVALUE=$(expr ${STASKID} - 1) #Split start at 0 and Task count start at 1. Substract 1 to Task count to match split name
DIGITID="0000${REALVALUE}"
DIGITID="${DIGITID: -4}"

echo "python ${SCRIPTDIR}/MakeAssignation.py -v ${VARNAME} -q ${SEQ_NUMBER} -1 ${PID}_${VARNAME}_R1.fastq -2 ${PID}_${VARNAME}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}"
python ${SCRIPTDIR}/MakeAssignation.py -v ${VARNAME} -q ${SEQ_NUMBER} -1 ${PID}_${VARNAME}_R1.fastq -2 ${PID}_${VARNAME}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}
