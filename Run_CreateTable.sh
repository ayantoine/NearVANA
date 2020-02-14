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
echo "python ${SDIR}/CreateTable.py -j ${nb_jobs} -p ${PID} -d ${DATA} -l ${VIRMINLEN}"
python ${SDIR}/CreateTable.py -j ${nb_jobs} -p ${PID} -d ${DATA} -l ${VIRMINLEN}
echo "------ /Create table ------"

echo "------ Xlsx conversion ------"
cat ${SDIR}/Tab2Xls.pl | wc -l
perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_Diamond_results.tab ${PID}_Diamond_results.xlsx $((${#PID}+7))
echo "------ /Xlsx conversion ------"

touch ${PID}.Table.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Table: "$delta > Time11.txt
