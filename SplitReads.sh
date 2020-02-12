#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
source $DATA

FASTQ=$2
PAIR=$3
VARNAME=$4
VAR_DODE="${VARNAME}[3]"

echo "------ Get Subsample list ------"
declare -a SAMPLE_LIST
while read c1 leftovers; do
	SAMPLE_LIST+=(${VARNAME}${c1})
done < ${!VAR_DODE}
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Subsample list ------"

SAMPLE=${SAMPLE_LIST[${STASKID}-1]}

echo "------ Split fastq by sample ------"
python ${SDIR}/SplitReads.py -f ${FASTQ} -r ${PID}_${VARNAME}_Hyper_Identified.tsv -s ${SAMPLE} -i ${PAIR}
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

