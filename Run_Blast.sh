#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Configure job array ------"
nb_seq=$(grep -c "^>" ${PID}_All.fa)
nb_task=${SMAXSIMJOB}
v_nb_seq_split=1000
let $[ nb_seq_split=$nb_seq/$nb_task ]
if [ $nb_seq_split == 0 ]
	then
	nb_seq_split=1
elif [ $nb_seq_split -gt $v_nb_seq_split ]
	then
	nb_seq_split=$v_nb_seq_split
fi

mkdir ${PID}_ToBlast
awk -v nb=$(($nb_seq_split*2)) 'NR%nb==1{file=FILENAME"."++i;}{print > file}' ${PID}_All.fa
mv ${PID}_All.fa.* ${PID}_ToBlast
nb_jobs=$(ls ${PID}_ToBlast | wc -l)

echo "Number of initial sequences: "$nb_seq
echo "Number of final sequences: "$(cat ${PID}"_ToBlast"/* | grep -c "^>")
echo "Number of sequences per job: "$nb_seq_split
echo "Number of generated jobs: "$nb_jobs
echo "Number of tasks: "$nb_task
echo "------ /Configure job array ------"


echo "------ Launch Blast by task ------"
echo "$SCALL $SPARAM $SRENAME ${PID}_N_Blast -e Run_BlastN.e -o Run_BlastN.o ${SDIR}/Run_BlastTask.sh $ARG N"
$SCALL $SPARAM $SRENAME ${PID}_N_Blast -e Run_BlastN.e -o Run_BlastN.o ${SDIR}/Run_BlastTask.sh $ARG N
echo "$SCALL $SPARAM $SRENAME ${PID}_N_Blast -e Run_BlastN.e -o Run_BlastN.o ${SDIR}/Run_BlastTask.sh $ARG X"
$SCALL $SPARAM $SRENAME ${PID}_X_Blast -e Run_BlastX.e -o Run_BlastX.o ${SDIR}/Run_BlastTask.sh $ARG X
while [ ! -e ${PID}.BlastN.ok ]; do sleep 60 ; done
while [ ! -e ${PID}.BlastX.ok ]; do sleep 60 ; done
echo "------ /Launch Blast by task ------"

rm -r ${PID}_ToBlast

touch ${PID}.Blast.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Blast: "$delta > Time09.txt
