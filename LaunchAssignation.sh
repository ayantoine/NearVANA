#! /bin/bash

KMER_FILE=$1
WORKDIR=$2
SCRIPTDIR=$3
OKDIR=$4
CONFFILE=$5

source ${CONFFILE}

REALVALUE=$(expr ${STASKID} - 1) #Split start at 0 and Task count start at 1. Substract 1 to Task count to match split name
DIGITID="000${REALVALUE}"
DIGITID="${DIGITID: -3}"


python ${SCRIPTDIR}/MakeAssignation.py -1 1R_${DIGITID} -2 2R_${DIGITID} -k ${KMER_FILE} -d ${WORKDIR} -t ${OKDIR}/${DIGITID}_MakeAssignation.ok
