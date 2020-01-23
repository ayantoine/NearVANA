#! /bin/bash

set -e

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

FASTQ=$2
PAIR=$3

echo "------ Get Sample list ------"
declare -a SAMPLE_LIST
while read c1 leftovers; do
	SAMPLE_LIST+=($c1)
done < ${DODE}
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Sample list ------"

SAMPLE=${SAMPLE_LIST[${STASKID}-1]}

echo "------ Split fastq by sample ------"
python ${SDIR}/SplitReads.py -f ${FASTQ} -r ${PID}_Hyper_Identified.tsv -s ${SAMPLE} -p ${PID} -i ${PAIR}
echo "------ /Split fastq by sample ------"

echo "Creating ok file"
while [ ! -f SplitReads${PAIR}_Ok/${STASKID}.SplitReads.${PAIR}.ok ]; do
	touch SplitReads${PAIR}_Ok/${STASKID}.SplitReads.${PAIR}.ok
	sleep 10
	if [ ! -d SplitReads${PAIR}_Ok ]; then
		break
	fi
done

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time SplitReads: "$delta

