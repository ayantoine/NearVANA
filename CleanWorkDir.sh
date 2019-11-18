#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

TASKARRAY=("N" "X")

echo "------ Compress All.fa ------"
gzip -f ${PID}_All.fa > ${PID}_All.fa.gz
echo "------ /Compress All.fa ------"

echo "------ Merge and compress rejected sequence ------"
for TASK in ${TASKARRAY[@]}; do
	touch ${PID}_Blast${TASK}.rejected.fa
	for FILE in ${PID}_BlastN/*.rejected ; do
		cat ${FILE} >> ${PID}_Blast${TASK}.rejected.fa
		sleep 1
		rm ${FILE}
	done
done
gzip -f ${PID}_Blast${TASK}.rejected.fa > ${PID}_Blast${TASK}.rejected.fa.gz
echo "------ /Merge and compress rejected sequence ------"

echo "------ Merge and compress keeped sequence ------"
for TASK in ${TASKARRAY[@]}; do
	touch ${PID}_Blast${TASK}.keeped.fa
	for FILE in ${PID}_BlastN/*.keeped ; do
		cat ${FILE} >> ${PID}_Blast${TASK}.keeped.fa
		sleep 1
		rm ${FILE}
	done
done
gzip -f ${PID}_Blast${TASK}.keeped.fa > ${PID}_Blast${TASK}.keeped.fa.gz
echo "------ /Merge and compress keeped sequence ------"

for TASK in ${TASKARRAY[@]}; do
	rm -r ${PID}_Blast${TASK}
done

touch ${PID}.Clean.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Assembly: "$delta > Time012.txt
