#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

nb_jobs=$(ls ${PID}_BlastX/*_2.tab | wc -l)

echo "------ Create short-list ------"
cut -f2,4 ${PID}_All.SPAdes_reverseAssembly.tsv | sort -u > ${PID}_All.SPAdes.contigs2sample.tsv
cut -f2,4 ${PID}_All.FLASH_reverseAssembly.tsv | sort -u > ${PID}_All.FLASH.contigs2sample.tsv
echo "------ /Create short-list ------"

echo "------ Create table ------"
	echo "$SCALL $SPARAM $SRENAME ${PID}_XTable -e Creation_TableX.e -o Creation_TableX.o ${SDIR}/CreateTable.sh $ARG X ${nb_jobs}"
	$SCALL $SPARAM $SRENAME ${PID}_XTable -e Creation_TableX.e -o Creation_TableX.o ${SDIR}/CreateTable.sh $ARG X ${nb_jobs}
	echo "$SCALL $SPARAM $SRENAME ${PID}_NTable -e Creation_TableN.e -o Creation_TableN.o ${SDIR}/CreateTable.sh $ARG N ${nb_jobs}"
	$SCALL $SPARAM $SRENAME ${PID}_NTable -e Creation_TableN.e -o Creation_TableN.o ${SDIR}/CreateTable.sh $ARG N ${nb_jobs}
	echo "$SCALL $SPARAM $SRENAME ${PID}_NTable -e Creation_TableAll.e -o Creation_TableAll.o ${SDIR}/CreateTableFusion.sh $ARG ${nb_jobs}"
	$SCALL $SPARAM $SRENAME ${PID}_AllTable -e Creation_TableAll.e -o Creation_TableAll.o ${SDIR}/CreateTableFusion.sh $ARG ${nb_jobs}
	while [ ! -e ${PID}.creationN.ok ]; do sleep 60 ; done
	while [ ! -e ${PID}.creationX.ok ]; do sleep 60 ; done
	while [ ! -e ${PID}.creationAll.ok ]; do sleep 60 ; done
	rm ${PID}.creationN.ok ${PID}.creationX.ok ${PID}.creationAll.ok
echo "------ /Create table ------"

echo "------ Xlsx conversion ------"
cat ${SDIR}/Tab2Xls.pl | wc -l
perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_BlastN_result.tab ${PID}_BlastN_result.xlsx $((${#PID}+7))
perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_BlastX_result.tab ${PID}_BlastX_result.xlsx $((${#PID}+7))
perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_BlastAll_result.tab ${PID}_BlastAll_result.xlsx $((${#PID}+7))
echo "------ /Xlsx conversion ------"

touch ${PID}.Table.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Table: "$delta > Time11.txt
