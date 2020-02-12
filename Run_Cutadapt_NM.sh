#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
source $DATA

echo "------ Get adaptors ------"
adapter1=$(cut -f1 ${ADAP})
adapter2=$(cut -f2 ${ADAP})
echo "Adapter1: "$adapter1
echo "Adapter2: "$adapter2
echo "------ /Get adaptors ------"

echo "------ Trim adaptors ------"
for VARNAME in "${PLATE[@]}"; do
	R1=${PID}_${VARNAME}_R1.fastq
	R2=${PID}_${VARNAME}_R2.fastq
	
	cutadapt --core=${MULTICPU} -a $adapter1 -A $adapter2 -q 30 -O $((${#adapter1}*85/100)) -m 15 -j 0 -o ${R1}.trim -p ${R2}.trim ${R1} ${R2}
	
	rm ${R1} ${R2}
done
echo "------ /Trim adaptors ------"

touch ${PID}.CutAdapt.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Trimming: "$delta > Time04.txt
