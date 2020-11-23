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

FASTQ=$2
PAIR=$3
VARNAME=$4
VAR_DODE="${VARNAME}[3]"

USE_KEEPUNASSIGNED="$(boolean "${UNASSIGNED}")"

echo "------ Get Subsample list ------"
declare -a SAMPLE_LIST
while read c1 leftovers; do
	SAMPLE_LIST+=(${VARNAME}${c1})
done < ${!VAR_DODE}
if [ "$USE_KEEPUNASSIGNED" = true ] ; then
	SAMPLE_LIST+=("UnassignedReads")
fi
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Subsample list ------"

SAMPLE=${SAMPLE_LIST[${STASKID}-1]}

echo "------ Split fastq by sample ------"
if  [ "$USE_KEEPUNASSIGNED" = true ] ; then
	echo "${SDIR}/SplitReads.py -u true -f ${FASTQ} -r ${PID}_${VARNAME}_Demultiplexing_Global.tsv -s ${SAMPLE} -i ${PAIR} -o ${SAMPLE}_${PID}_R${PAIR}.fastq.split"
	python ${SDIR}/SplitReads.py -u true -f ${FASTQ} -r ${PID}_${VARNAME}_Demultiplexing_Global.tsv -s ${SAMPLE} -i ${PAIR} -o ${SAMPLE}_${PID}_R${PAIR}.fastq.split
else
	echo "${SDIR}/SplitReads.py -u true -f ${FASTQ} -r ${PID}_${VARNAME}_Demultiplexing_Global.tsv -s ${SAMPLE} -i ${PAIR} -o ${SAMPLE}_${PID}_R${PAIR}.fastq.split"
	python ${SDIR}/SplitReads.py -f ${FASTQ} -r ${PID}_${VARNAME}_Hyper_Identified.tsv -s ${SAMPLE} -i ${PAIR} -o ${SAMPLE}_${PID}_R${PAIR}.fastq.split
fi
echo "------ /Split fastq by sample ------"

echo "Creating ok file"
while [ ! -f SplitReads${VARNAME}-${PAIR}_Ok/${STASKID}.SplitReads.${VARNAME}-${PAIR}.ok ]; do
	touch SplitReads${VARNAME}-${PAIR}_Ok/${STASKID}.SplitReads.${VARNAME}-${PAIR}.ok
	sleep 10
	if [ ! -d SplitReads${VARNAME}-${PAIR}_Ok ]; then
		break
	fi
done

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time SplitReads: "$delta

