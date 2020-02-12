#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
source $DATA

echo "------ Get Subsample list ------"
declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	VAR_SAMPLE_FILE="${VARNAME}[3]"
	#echo "${!VAR_SAMPLE_FILE}"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Subsample list ------"

SAMPLE=${SAMPLE_LIST[${STASKID}-1]}

echo "------ Retrieve Pair in sample ------"
echo "python ${SDIR}/RetrievePair.py -i ${SAMPLE}/${SAMPLE}_${PID}_R1.fastq.split.trim -p ${SAMPLE}/${SAMPLE}_${PID}_R2.fastq.split.trim"
python ${SDIR}/RetrievePair.py -i ${SAMPLE}/${SAMPLE}_${PID}_R1.fastq.split.trim -p ${SAMPLE}/${SAMPLE}_${PID}_R2.fastq.split.trim
rm ${SAMPLE}/${SAMPLE}_${PID}_R1.fastq.split.trim ${SAMPLE}/${SAMPLE}_${PID}_R2.fastq.split.trim
echo "------ /Retrieve Pair in sample ------"

touch RetrievePair_Ok/${STASKID}.RetrievePair.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrivePair: "$delta

