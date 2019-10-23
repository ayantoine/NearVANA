#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

FASTQ=$2
PAIR=$3

echo "------ Launch SplitReads array ------"
nb_jobs=$(cut -f1 ${DODE} | wc -l)
if [ ! -d "SplitReads${PAIR}_Ok" ] ; then mkdir "SplitReads${PAIR}_Ok" ; fi
if [ ! -d ${PID}"_log_SplitReads${PAIR}" ] ; then mkdir ${PID}"_log_SplitReads${PAIR}" ; fi
echo "$SCALL $SPARAM $SRENAME ${PID}_${PAIR}-SplitReads ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_SplitReads${PAIR}"/${PID}_SplitReads${PAIR}.e${SPSEUDOTASKID} -o ${PID}"_log_SplitReads${PAIR}"/${PID}_SplitReads${PAIR}.o${SPSEUDOTASKID} ${SDIR}/SplitReads.sh ${ARG} ${FASTQ} ${PAIR}"
$SCALL $SPARAM $SRENAME ${PID}_${PAIR}-SplitReads ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_SplitReads${PAIR}"/${PID}_SplitReads${PAIR}.e${SPSEUDOTASKID} -o ${PID}"_log_SplitReads${PAIR}"/${PID}_SplitReads${PAIR}.o${SPSEUDOTASKID} ${SDIR}/SplitReads.sh ${ARG} ${FASTQ} ${PAIR}
while true ; do
	if [ $(ls SplitReads${PAIR}_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls SplitReads${PAIR}_Ok/*_SplitReads${PAIR}.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		then
		rm -r SplitReads${PAIR}_Ok
		break
	fi
	sleep 60
done
echo "------ /Launch SplitReads array ------"

> ${FASTQ}.split.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Demultiplexing: "$delta > Time03-${PAIR}.txt

