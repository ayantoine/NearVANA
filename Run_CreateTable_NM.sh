#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

nb_jobs=$(ls ${PID}_Diamond/*_2.tab | wc -l)

echo "------ Create short-list ------"
cut -f2,4 ${PID}_All.Megahit_reverseAssembly.tsv | sort -u > ${PID}_All.Megahit.contigs2sample.tsv
cut -f2,4 ${PID}_All.FLASH_reverseAssembly.tsv | sort -u > ${PID}_All.FLASH.contigs2sample.tsv
echo "------ /Create short-list ------"

echo "------ Create table ------"
	echo "$SCALL $SPARAM $SRENAME ${PID}_Table -e Creation_Table.e -o Creation_Table.o ${SDIR}/CreateTable_NM.sh $ARG ${nb_jobs}"
	$SCALL $SPARAM $SRENAME ${PID}_Table -e Creation_Table.e -o Creation_Table.o ${SDIR}/CreateTable_NM.sh $ARG ${nb_jobs}
	while [ ! -e ${PID}.creationTable.ok ]; do sleep 60 ; done
	rm ${PID}.creationTable.ok
echo "------ /Create table ------"

echo "------ Xlsx conversion ------"
cat ${SDIR}/Tab2Xls.pl | wc -l
perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_Diamond_results.tab ${PID}_Diamond_results.xlsx $((${#PID}+7))
echo "------ /Xlsx conversion ------"

touch ${PID}.Table.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Table: "$delta > Time11.txt
