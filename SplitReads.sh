#! /bin/bash

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
python ${SDIR}/SplitReads.py -f ${FASTQ} -r ${PID}_Hyper_Identified.tab -s ${SAMPLE} -p ${PID} -i ${PAIR}
echo "------ /Split fastq by sample ------"

touch SplitReads${PAIR}_Ok/${STASKID}.SplitReads.${PAIR}.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Demultiplexing: "$delta

