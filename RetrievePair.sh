#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Get Sample list ------"
declare -a SAMPLE_LIST
while read c1 leftovers; do
	SAMPLE_LIST+=($c1)
done < ${DODE}
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Sample list ------"

SAMPLE=${SAMPLE_LIST[${STASKID}-1]}

echo "------ Retrieve Pair in sample ------"
python ${SDIR}/RetrievePair.py -i ${STASKID}/${STASKID}_${PID}_R1.fastq.split.trim -p ${STASKID}/${STASKID}_${PID}_R2.fastq.split.trim
echo "------ /Retrieve Pair in sample ------"

touch RetrivePair_Ok/${STASKID}.RetrivePair.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrivePair: "$delta

