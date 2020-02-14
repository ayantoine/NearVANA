#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

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

mkdir ${PID}_ToDiamond

echo "python ${SDIR}/SplitFasta.py -f ${PID}_ToDiamond -i ${PID}_All.fa c ${CHUNCK}"
python ${SDIR}/SplitFasta.py -f ${PID}_ToDiamond -i ${PID}_All.fa -c ${CHUNCK}

nb_jobs=$(ls ${PID}_ToDiamond | wc -l)

echo "Number of initial sequences: "$nb_seq
echo "Number of final sequences: "$(cat ${PID}"_ToDiamond"/* | grep -c -f Target.txt)
echo "Number of sequences per job: "$CHUNCK
echo "Number of generated jobs: "$nb_jobs
echo "Number of simultaneous tasks: "$nb_task
echo "------ /Configure job array ------"
rm Target.txt

echo "------ Launch Diamond ------"
echo "$SCALL $SPARAM $SRENAME ${PID}_Diamond -e Run_Diamond.e -o Run_Diamond.o ${SDIR}/Diamond.sh $ARG"
$SCALL $SPARAM $SRENAME ${PID}_DiamondTask -e Run_DiamondTask.e -o Run_DiamondTask.o ${SDIR}/Run_DiamondTask.sh $ARG
while [ ! -e ${PID}.DiamondTask.ok ]; do sleep 60 ; done
echo "------ Launch /Diamond ------"

rm -r ${PID}_ToDiamond

echo "------ Compress All.fa ------"
gzip -f ${PID}_All.fa > ${PID}_All.fa.gz
echo "------ /Compress All.fa ------"

touch ${PID}.Diamond.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Diamond: "$delta > Time09.txt
