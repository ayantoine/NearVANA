#!/bin/bash

ARG=$1
source $ARG
source $CONF

REALVALUE=$((${STASKID} - 1)) #SGE count from 1, linux cut count from 0

DIGITID="000${REALVALUE}"
DIGITID="${DIGITID: -3}"

adapter1=$(cut -f1 ${$ADAP})
adapter2=$(cut -f2 ${$ADAP})
echo
echo "----------------------------------------"
echo "Adapter1: "$adapter1
echo "Adapter2: "$adapter2
echo "----------------------------------------"
cutadapt -a $adapter1 -A $adapter2 -q 30 -O $((${#adapter1}*85/100)) -m 15 -j 0 -o ${PID}_Cleaning/${DIGITID}_R1.Trim1.fastq -p ${PID}_Cleaning/${DIGITID}_R2.Trim1.fastq ${PID}_Cleaning/1Clean_${DIGITID} ${PID}_Cleaning/2Clean_${DIGITID}

touch CutAdapt_Ok/${DIGITID}_CutAdapat.ok
