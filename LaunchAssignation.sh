#! /bin/bash

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

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

USE_PAIREND="$(boolean "${PAIREND}")"

REALVALUE=$(expr ${STASKID} - 1) #Split start at 0 and Task count start at 1. Substract 1 to Task count to match split name
DIGITID="0000${REALVALUE}"
DIGITID="${DIGITID: -4}"

if [ "$USE_PAIREND" = true ] ; then
	echo "python ${SCRIPTDIR}/MakeAssignation.py -p 1 -v ${VARNAME} -q ${SEQ_NUMBER} -1 ${PID}_${VARNAME}_R1.fastq -2 ${PID}_${VARNAME}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}"
	python ${SCRIPTDIR}/MakeAssignation.py -p True -v ${VARNAME} -q ${SEQ_NUMBER} -1 ${PID}_${VARNAME}_R1.fastq -2 ${PID}_${VARNAME}_R2.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}
else
	echo "python ${SCRIPTDIR}/MakeAssignation.py -p 0 -v ${VARNAME} -q ${SEQ_NUMBER} -1 ${PID}_${VARNAME}_R1.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}"
	python ${SCRIPTDIR}/MakeAssignation.py -p False -v ${VARNAME} -q ${SEQ_NUMBER} -1 ${PID}_${VARNAME}_R1.fastq -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok -i ${DIGITID} -c ${CONFFILE}
fi
