#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

FASTQ=$2
TASK=$3
PAIR=$4
OK_FOLDER=$5

echo "------ Get Sample list ------"
declare -a SAMPLE_LIST
while read c1 leftovers; do
	SAMPLE_LIST+=($c1)
done < ${DODE}
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Sample list ------"

SAMPLE=${SAMPLE_LIST[${TASK}-1]}

echo "------ Split fastq by sample ------"
python ${SDIR}/SplitRead.py -f ${FASTQ} -r ${PID}_Hyper_Identified.tab -s ${SAMPLE} -p ${PID} -i ${PAIR}
echo "------ /Split fastq by sample ------"

> ${OK_FOLDER}/${TASK}.${FASTQ}.split.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Demultiplexing: "$delta > Time02.txt

