#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Launch RetrivePair array ------"
nb_jobs=$(cut -f1 ${DODE} | wc -l)
if [ ! -d "RetrivePair_Ok" ] ; then mkdir "RetrivePair_Ok" ; fi
if [ ! -d ${PID}"_log_RetrivePair" ] ; then mkdir ${PID}"_log_RetrivePair" ; fi
echo "$SCALL $SPARAM $SRENAME ${PID}_${PAIR}-RetrievePair ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_RetrievePair"/${PID}_RetrievePair.e${SPSEUDOTASKID} -o ${PID}"_log_RetrievePair"/${PID}_RetrievePair.o${SPSEUDOTASKID} ${SDIR}/RetrievePair.sh ${ARG}"
$SCALL $SPARAM $SRENAME ${PID}_${PAIR}-RetrievePair ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_RetrievePair"/${PID}_RetrievePair.e${SPSEUDOTASKID} -o ${PID}"_log_RetrievePair"/${PID}_RetrievePair.o${SPSEUDOTASKID} ${SDIR}/RetrievePair.sh ${ARG}
while true ; do
	if [ $(ls RetrivePair_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls RetrivePair_Ok/*.RetrivePair.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		then
		rm -r RetrivePair_Ok
		break
	fi
	sleep 60
done
echo "------ /Launch RetrivePair array ------"

touch ${PID}.fastq.deinterlaced.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrievePair: "$delta > Time05.txt

