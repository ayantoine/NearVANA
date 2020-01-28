#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

#N or X
TASK=$2

nb_jobs=$(ls ${PID}_ToBlast | wc -l)

if [ ${TASK} == N ]; then
    TAXO=${NUCTAXO}
    TIMERID="1"
elif [ ${TASK} == X ]; then
    TAXO=${PROTAXO}
    TIMERID="2"
fi

if [ ! -d "Blast${TASK}_Ok" ] ; then mkdir "Blast${TASK}_Ok" ; fi
if [ ! -d ${PID}"_Blast${TASK}" ] ; then mkdir ${PID}"_Blast${TASK}" ; fi
if [ ! -d ${PID}"_log_Blast${TASK}" ] ; then mkdir ${PID}"_log_Blast${TASK}" ; fi

if [ ! -f ${PID}.Blast${TASK}.ok ] ; then
    echo "$SCALL $SPARAM $SRENAME ${PID}_${TASK}-Blast${TASK} ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Blast${TASK}"/${PID}_Blast${TASK}.e${SPSEUDOTASKID} -o ${PID}"_log_Blast${TASK}"/${PID}_Blast${TASK}.o${SPSEUDOTASKID} ${SDIR}/CutBlastJob.sh $ARG ${TASK}"
    $SCALL $SPARAM $SRENAME ${PID}_${TASK}-Blast${TASK} ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Blast${TASK}"/${PID}_Blast${TASK}.e${SPSEUDOTASKID} -o ${PID}"_log_Blast${TASK}"/${PID}_Blast${TASK}.o${SPSEUDOTASKID} ${SDIR}/BlastJob.sh $ARG ${TASK} # Input: ${PID}_All_sequences.fa ${BK} ${PID} ${SDIR} ${PROTAXO}
    #while true ; do
	    #if [ $(ls Blast${TASK}_Ok/ | wc -l) -eq 0 ]
		    #then
		    #nbr_ok=0
	    #else
		    #nbr_ok=$(ls Blast${TASK}_Ok/*_Blast${TASK}.ok | wc -l)
	    #fi
	    #if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		    #then
		    #rm -r Blast${TASK}_Ok
		    #break
	    #fi
	    #sleep 60
    #done
    rm -r Blast${TASK}_Ok
else
    echo "${PID}.Blast${TASK}.ok already existing, do nothing..."
fi

touch ${PID}.Blast${TASK}.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Blast${TASK}: "$delta > Time09-${TIMERID}.txt
