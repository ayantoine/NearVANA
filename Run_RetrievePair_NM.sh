#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
source $DATA

echo "------ Retrieve Pair in sample ------"
for VARNAME in "${PLATE[@]}"; do
	R1=${PID}_${VARNAME}_R1.fastq.trim
	R2=${PID}_${VARNAME}_R2.fastq.trim
	
	echo "python ${SDIR}/RetrievePair.py -i ${R1} -p ${R2}"
	python ${SDIR}/RetrievePair.py -i ${R1} -p ${R2}
	
	rm ${R1} ${R2}
done
echo "------ /Retrieve Pair in sample ------"

touch ${PID}.Deinterlacing.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrievePair: "$delta > Time05.txt

