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

mkdir ${PID}_ToBlast

echo "python ${SDIR}/SplitFasta.py -f ${PID}_ToBlast -i ${PID}_All.fa"
python ${SDIR}/SplitFasta.py -f ${PID}_ToBlast -i ${PID}_All.fa

nb_jobs=$(ls ${PID}_ToBlast | wc -l)

echo "Number of initial sequences: "$nb_seq
echo "Number of final sequences: "$(cat ${PID}"_ToBlast"/* | grep -c -f Target.txt)
echo "Number of sequences per job: "$(cat ${PID}"_ToBlast"/${PID}"_All.fa.1" | grep -c -f Target.txt)
echo "Number of generated jobs: "$nb_jobs
echo "Number of simultaneous tasks: "$nb_task
echo "------ /Configure job array ------"
rm Target.txt


echo "------ Launch Blast by task ------"
echo "$SCALL $SPARAM $SRENAME ${PID}_N_Blast -e Run_BlastN.e -o Run_BlastN.o ${SDIR}/Run_BlastTask.sh $ARG N"
$SCALL $SPARAM $SRENAME ${PID}_N_Blast -e Run_BlastN.e -o Run_BlastN.o ${SDIR}/Run_BlastTask.sh $ARG N
echo "$SCALL $SPARAM $SRENAME ${PID}_X_Blast -e Run_BlastX.e -o Run_BlastX.o ${SDIR}/Run_BlastTask.sh $ARG X"
$SCALL $SPARAM $SRENAME ${PID}_X_Blast -e Run_BlastX.e -o Run_BlastX.o ${SDIR}/Run_BlastTask.sh $ARG X
while [ ! -e ${PID}.BlastN.ok ]; do sleep 60 ; done
while [ ! -e ${PID}.BlastX.ok ]; do sleep 60 ; done
echo "------ /Launch Blast by task ------"

rm -r ${PID}_ToBlast

echo "------ Compress All.fa ------"
gzip -f ${PID}_All.fa > ${PID}_All.fa.gz
echo "------ /Compress All.fa ------"

touch ${PID}.Blast.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Blast: "$delta > Time09.txt
