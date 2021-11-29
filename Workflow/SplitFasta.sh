#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

echo "------ Configure job array ------"
SEQMARKER="^>"
echo ${SEQMARKER} > Target.txt
nb_seq=$(grep -c -f Target.txt ${PID}_All.fa)
nb_task=${SMAXSIMJOB}

echo "$SMAXARRAYSIZE"
if [[ $SMAXARRAYSIZE -eq 0 ]]; then
	CHUNCK=1000
else
	CHUNCK=$((nb_seq/SMAXARRAYSIZE+1))
fi

mkdir ${PID}_ToBlast

echo "python ${SDIR}/SplitFasta.py -f ${PID}_ToBlast -i ${PID}_All.fa c ${CHUNCK}"
python ${SDIR}/SplitFasta.py -f ${PID}_ToBlast -i ${PID}_All.fa -c ${CHUNCK}

nb_jobs=$(ls ${PID}_ToBlast | wc -l)

echo "Number of initial sequences: "$nb_seq
echo "Number of final sequences: "$(cat ${PID}"_ToBlast"/* | grep -c -f Target.txt)
echo "Number of sequences per job: "$CHUNCK
echo "Number of generated jobs: "$nb_jobs
echo "Number of simultaneous tasks: "$nb_task
echo "------ /Configure job array ------"
rm Target.txt

echo "------ Compress All.fa ------"
gzip -f ${PID}_All.fa > ${PID}_All.fa.gz
echo "------ /Compress All.fa ------"

touch ${PID}.SplitFasta.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time SplitData: "$delta > Time08.txt

