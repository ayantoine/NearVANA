#! /bin/bash

datetime1=$(date +%s)

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

ARG=$1
source $ARG
source $CONF
source $DATA
SDIR=${GITDIR}/Workflow

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

echo "------ Retrieve Pair in sample ------"
for VARNAME in "${PLATE[@]}"; do

	if [ "$USE_PAIREND" = true ] ; then
		R1=${PID}_${VARNAME}_R1.fastq.trim
		R2=${PID}_${VARNAME}_R2.fastq.trim
	else
		R1=${PID}_${VARNAME}_R1.fastq.trim
		R2=${PID}_${VARNAME}_R2.fastq.trim
		touch ${PID}_${VARNAME}_R2.fastq.trim
	fi
	echo "python ${SDIR}/RetrievePair.py -i ${R1} -p ${R2}"
	python ${SDIR}/RetrievePair.py -i ${R1} -p ${R2}
	
	rm ${R1} ${R2}
done
echo "------ /Retrieve Pair in sample ------"

touch ${PID}.Deinterlacing.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrievePair: "$delta > Time05.txt

