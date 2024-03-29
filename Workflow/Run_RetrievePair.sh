#! /bin/bash

datetime1=$(date +%s)

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

ARG=$1
source $ARG
source $CONF
source $DATA
SDIR=${GITDIR}/Workflow

USE_PAIREND="$(boolean "${PAIREND}")"
USE_METADATA="$(boolean "${METADATA}")"
USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"
USE_KEEPUNASSIGNED="$(boolean "${UNASSIGNED}")"

NB_ITEM=1
ID_R1=0
if [ "$USE_PAIREND" = true ] ; then
	ID_R2=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi
if [ "$USE_MULTIPLEX" = true ] ; then
	ID_DODE=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi
if [ "$USE_METADATA" = true ] ; then
	ID_META=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi

echo "------ Get Subsample list ------"
declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	VAR_SAMPLE_FILE="${VARNAME}[$ID_DODE]"
	#echo "${!VAR_SAMPLE_FILE}"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done
if [ "$USE_KEEPUNASSIGNED" = true ] ; then
	SAMPLE_LIST+=("UnassignedReads")
fi
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Subsample list ------"

echo "------ Launch RetrievePair array ------"
nb_jobs=${#SAMPLE_LIST[@]} #$(cut -f1 ${DODE} | wc -l)
if [ ! -d "RetrievePair_Ok" ] ; then mkdir "RetrievePair_Ok" ; fi
if [ ! -d ${PID}"_log_RetrievePair" ] ; then mkdir ${PID}"_log_RetrievePair" ; fi
echo "$SCALL $SPARAM_HEAVY $SRENAME ${PID}_RetrievePair ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_RetrievePair"/${PID}_RetrievePair.e${SPSEUDOTASKID} -o ${PID}"_log_RetrievePair"/${PID}_RetrievePair.o${SPSEUDOTASKID} ${SDIR}/RetrievePair.sh ${ARG}"
$SCALL $SPARAM_HEAVY $SRENAME ${PID}_RetrievePair ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_RetrievePair"/${PID}_RetrievePair.e${SPSEUDOTASKID} -o ${PID}"_log_RetrievePair"/${PID}_RetrievePair.o${SPSEUDOTASKID} ${SDIR}/RetrievePair.sh ${ARG}
while true ; do
	if [ $(ls RetrievePair_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls RetrievePair_Ok/*.RetrievePair.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		then
		rm -r RetrievePair_Ok
		break
	fi
	sleep 60
done
echo "------ /Launch RetrievePair array ------"

touch ${PID}.Deinterlacing.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time RetrievePair: "$delta > Time05.txt

