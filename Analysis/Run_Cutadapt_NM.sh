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


echo "------ Get adaptors ------"
adapter1=$(cut -f1 ${ADAP})
adapter2=$(cut -f2 ${ADAP})
echo "Adapter1: "$adapter1
echo "Adapter2: "$adapter2
echo "------ /Get adaptors ------"

echo "------ Trim adaptors ------"
for VARNAME in "${PLATE[@]}"; do
	R1=${PID}_${VARNAME}_R1.fastq
	if [ "$USE_PAIREND" = true ] ; then
		R2=${PID}_${VARNAME}_R2.fastq
		cutadapt --core=${MULTICPU} -a $adapter1 -A $adapter2 -q 30 -O $((${#adapter1}*85/100)) -m 15 -j 0 -o ${R1}.trim -p ${R2}.trim ${R1} ${R2}
		rm ${R1} ${R2}
	else
		cutadapt --core=${MULTICPU} -a $adapter1 -q 30 -O $((${#adapter1}*85/100)) -m 15 -j 0 -o ${R1}.trim ${R1}
		rm ${R1}
	fi	
done
echo "------ /Trim adaptors ------"

touch ${PID}.CutAdapt.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Trimming: "$delta > Time04.txt
