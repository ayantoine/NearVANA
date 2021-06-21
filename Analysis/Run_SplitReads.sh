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
SDIR=${GITDIR}/Analysis

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

FASTQ=$2
PAIR=$3
VARNAME=$4
VAR_DODE="${VARNAME}[$ID_DODE]"

echo "------ Launch SplitReads array ------"
nb_jobs=$(cut -f1 ${!VAR_DODE} | wc -l)
if [ "$USE_KEEPUNASSIGNED" = true ] ; then
	let nb_jobs++
fi

if [ ! -d "SplitReads${VARNAME}-${PAIR}_Ok" ] ; then mkdir "SplitReads${VARNAME}-${PAIR}_Ok" ; fi
if [ ! -d ${PID}"_log_SplitReads${VARNAME}-${PAIR}" ] ; then mkdir ${PID}"_log_SplitReads${VARNAME}-${PAIR}" ; fi
echo "$SCALL $SPARAM $SRENAME ${PID}_${PAIR}-SplitReads ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_SplitReads${VARNAME}-${PAIR}"/${PID}_SplitReads${VARNAME}-${PAIR}.e${SPSEUDOTASKID} -o ${PID}"_log_SplitReads${VARNAME}-${PAIR}"/${PID}_SplitReads${VARNAME}-${PAIR}.o${SPSEUDOTASKID} ${SDIR}/SplitReads.sh ${ARG} ${FASTQ} ${PAIR} ${VARNAME}"
 $SCALL $SPARAM $SRENAME ${PID}_${PAIR}-SplitReads ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_SplitReads${VARNAME}-${PAIR}"/${PID}_SplitReads${VARNAME}-${PAIR}.e${SPSEUDOTASKID} -o ${PID}"_log_SplitReads${VARNAME}-${PAIR}"/${PID}_SplitReads${VARNAME}-${PAIR}.o${SPSEUDOTASKID} ${SDIR}/SplitReads.sh ${ARG} ${FASTQ} ${PAIR} ${VARNAME}
while true ; do
	if [ $(ls SplitReads${VARNAME}-${PAIR}_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls SplitReads${VARNAME}-${PAIR}_Ok/*.SplitReads.${VARNAME}-${PAIR}.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		then
		rm -r SplitReads${VARNAME}-${PAIR}_Ok
		break
	fi
	sleep 60
done
echo "------ /Launch SplitReads array ------"

touch ${FASTQ}.split.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time SplitReads & Trim linkers: "$delta > Time03-${VARNAME}-${PAIR}.txt

