#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $DATA

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

USE_PAIREND="$(boolean "${PAIREND}")"
USE_METADATA="$(boolean "${METADATA}")"
USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"
USE_KEEPUNASSIGNED="$(boolean "${UNASSIGNED}")"

NB_ITEM=1
ID_R1=0
if [ "$USE_PAIREND" = true ] ; then
	ID_R2=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi
if [ "$USE_MULTIPLEX" = true ] ; then
	ID_DODE=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi
if [ "$USE_METADATA" = true ] ; then
	ID_META=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi

for VARNAME in "${PLATE[@]}"; do
	VAR_R1_FILE="${VARNAME}[$ID_R1]"
	zcat ${!VAR_R1_FILE} > ${PID}_${VARNAME}_R1.fastq
	
	if [ "$USE_PAIREND" = true ] ; then
		VAR_R2_FILE="${VARNAME}[$ID_R2]"
		zcat ${!VAR_R2_FILE} > ${PID}_${VARNAME}_R2.fastq
	fi
done

touch ${PID}.extraction.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Gz_extraction: "$delta > Time01.txt
