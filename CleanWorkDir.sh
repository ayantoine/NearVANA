#! /bin/bash

set -e

err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

TASKARRAY=("N" "X")

echo "------ Merge and compress rejected sequence ------"
for TASK in ${TASKARRAY[@]}; do
	touch ${PID}_Blast${TASK}.rejected.fa
	for FILE in ${PID}_Blast${TASK}/*.rejected ; do
		cat ${FILE} >> ${PID}_Blast${TASK}.rejected.fa
		sleep 1
		rm ${FILE}
	done
	sleep 1
	gzip -f ${PID}_Blast${TASK}.rejected.fa > ${PID}_Blast${TASK}.rejected.fa.gz
done
echo "------ /Merge and compress rejected sequence ------"

echo "------ Merge and compress keeped sequence ------"
for TASK in ${TASKARRAY[@]}; do
	touch ${PID}_Blast${TASK}.keeped.fa
	for FILE in ${PID}_Blast${TASK}/*.keeped ; do
		cat ${FILE} >> ${PID}_Blast${TASK}.keeped.fa
		sleep 1
		rm ${FILE}
	done
	sleep 1
	#gzip -f ${PID}_Blast${TASK}.keeped.fa > ${PID}_Blast${TASK}.keeped.fa.gz
done
echo "------ /Merge and compress keeped sequence ------"

echo "------ Compute stat ------"
zcat ${PID}_All.fa.gz > ${PID}_All.fa
python ${SDIR}/ComputeStat.py -p ${PID}
rm ${PID}_All.fa
gzip -f ${PID}_BlastX.keeped.fa > ${PID}_BlastX.keeped.fa.gz
gzip -f ${PID}_BlastN.keeped.fa > ${PID}_BlastN.keeped.fa.gz
echo "------ Compute stat ------"

for TASK in ${TASKARRAY[@]}; do
	rm -r ${PID}_Blast${TASK}
done

touch ${PID}.Clean.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Cleaning: "$delta > Time12.txt
