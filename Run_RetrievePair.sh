#! /bin/bash

set -e

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Launch RetrievePair array ------"
nb_jobs=$(cut -f1 ${DODE} | wc -l)
if [ ! -d "RetrievePair_Ok" ] ; then mkdir "RetrievePair_Ok" ; fi
if [ ! -d ${PID}"_log_RetrievePair" ] ; then mkdir ${PID}"_log_RetrievePair" ; fi
echo "$SCALL $SPARAM $SRENAME ${PID}_RetrievePair ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_RetrievePair"/${PID}_RetrievePair.e${SPSEUDOTASKID} -o ${PID}"_log_RetrievePair"/${PID}_RetrievePair.o${SPSEUDOTASKID} ${SDIR}/RetrievePair.sh ${ARG}"
$SCALL $SPARAM $SRENAME ${PID}_RetrievePair ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_RetrievePair"/${PID}_RetrievePair.e${SPSEUDOTASKID} -o ${PID}"_log_RetrievePair"/${PID}_RetrievePair.o${SPSEUDOTASKID} ${SDIR}/RetrievePair.sh ${ARG}
#while true ; do
	#if [ $(ls RetrievePair_Ok/ | wc -l) -eq 0 ]
		#then
		#nbr_ok=0
	#else
		#nbr_ok=$(ls RetrievePair_Ok/*.RetrievePair.ok | wc -l)
	#fi
	#if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		#then
		#rm -r RetrievePair_Ok
		#break
	#fi
	#sleep 60
#done
rm -r RetrievePair_Ok
echo "------ /Launch RetrievePair array ------"

touch ${PID}.Deinterlacing.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrievePair: "$delta > Time05.txt

