#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

nb_jobs=$(ls ${PID}_ToBlast | wc -l)

if [ ! -d "Diamond_Ok" ] ; then mkdir "Diamond_Ok" ; fi
if [ ! -d ${PID}"_BlastD" ] ; then mkdir ${PID}"_BlastD" ; fi
if [ ! -d ${PID}"_log_Diamond" ] ; then mkdir ${PID}"_log_Diamond" ; fi

if [ ! -f ${PID}.BlastD.ok ] ; then
    echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Diamond ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Diamond"/${PID}_Diamond.e${SPSEUDOTASKID} -o ${PID}"_log_Diamond"/${PID}_Diamond.o${SPSEUDOTASKID} ${SDIR}/Diamond.sh $ARG"
    $SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Diamond ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Diamond"/${PID}_Diamond.e${SPSEUDOTASKID} -o ${PID}"_log_Diamond"/${PID}_Diamond.o${SPSEUDOTASKID} ${SDIR}/Diamond.sh $ARG
    while true ; do
	    if [ $(ls Diamond_Ok/ | wc -l) -eq 0 ]
		    then
		    nbr_ok=0
	    else
		    nbr_ok=$(ls Diamond_Ok/*_Diamond.ok | wc -l)
	    fi
	    if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		    then
		    rm -r Diamond_Ok
		    break
	    fi
	    sleep 60
    done
else
    echo "${PID}.Diamond.ok already existing, do nothing..."
fi

touch ${PID}.BlastD.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Diamond: "$delta > Time09-D.txt
