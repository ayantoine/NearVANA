#! /bin/bash

ARG=$1
#WARNING!! Source file is executed (Security, etc.)
source $ARG

echo "------ Check Input existence ------"
if [ ! -d $SDIR ] ; then
	echo "Directory SDIR ${SDIR} does not exists"
	exit 1
fi
LIST_FILE=($R1 $R2 $ADAP $DODE $META $SUBS $NUCACC $NUCDEF $PROACC $PRODEF $DBLINEAGE $CONF)
for i in "${LIST_FILE[@]}"; do
	if [ ! -f $i ]; then
		echo "File $i does not exists"
		exit 1
	fi
done

echo "> Args details"
List_NONFILE=(PID R1 R2 ADAP DODE META SUBS NUCACC NUCDEF PROACC PRODEF DBLINEAGE CONF)
for i in "${List_NONFILE[@]}"; do
	echo "$i: ${!i}"
done

echo
echo '> Conf details (/!\ beware bash interpretation for STASKID)'
source $CONF
LIST_PARAM=(SCALL SPARAM STASKARRAY SMAXTASK SRENAME SMAXSIMJOB STASKID SPSEUDOTASKID VIRNTDB ALLNTDB VIRPTDB ALLPTDB)
for i in "${LIST_PARAM[@]}"; do
	echo "$i: ${!i}"
done
echo "------ /Check Input existence ------"


echo "------ Extract .gz ------"
if [ ! -f ${PID}_R1.fastq ] || [ ! -f ${PID}_R1.fastq ]; then
	echo "bash ${SDIR}/Gz_extraction.sh $ARG"
	bash ${SDIR}/Gz_extraction.sh $ARG # Output: ${PID}_R1.fastq ${PID}_R2.fastq
else
	echo "${PID}_R1.fastq and ${PID}_R2.fastq already existing, pass"
fi
echo "------ /Extract .gz ------"

echo "------ Demultiplexing reads ------"
if [ ! -f ${PID}_Hyper_Identified.tab ] || [ ! -f ${PID}_Demultiplexing_Hyper.tab ] || [ ! -f ${PID}_Demultiplexing_Hyper_Distribution.tab ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing_Illumina_pe_V5.e -o Demultiplexing_Illumina_pe_V5.o ${SDIR}/Demultiplexing_Illumina_pe_V5.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing.e -o Demultiplexing.o ${SDIR}/Demultiplexing.sh $ARG # Input: ${PID}_R1.fastq ${PID}_R2.fastq $MID $PID $SDIR # Output: ${PID}_Demultiplexing.tab ${PID}_Demultiplexing_Distribution.tab ${PID}_Hyper_Identified.tab ${PID}_Hypo_1_Identified.tab ${PID}_Hypo_2_Identified.tab ${PID}_Ambiguous.tab ${PID}_Unidentified.tab
else
	echo "${PID}_Hyper_Identified.tab, ${PID}_Demultiplexing_Hyper.tab and ${PID}_Demultiplexing_Hyper_Distribution.tab already existing, pass"
fi
echo "------ Demultiplexing reads Launched------"

echo "------ Trim adapters ------"
if [ ! -f ${PID}_R1.Trimmed.fastq ] || [ ! -f ${PID}_R2.Trimmed.fastq ]; then
	echo $SCALL $SPARAM $SRENAME ${PID}_Cleaning -e Cleaning.e -o Cleaning.o ${SDIR}/Cleaning.sh $ARG
	$SCALL $SPARAM $SRENAME ${PID}_Cleaning -e Cleaning.e -o Cleaning.o ${SDIR}/Cleaning.sh $ARG # Input: ${PID}_R1.fastq ${PID}_R2.fastq $MID $PID Output: ${PID}_R1.Trimmed.fastq ${PID}_R2.Trimmed.fastq
else
	echo "${PID}_R1.Trimmed.fastq and ${PID}_R2.Trimmed.fastq  already existing, pass"
fi
while [ ! -e ${PID}_Cleaning.ok ]; do sleep 60 ; done
echo "------ /Trim adapters ------"



echo "------ Waiting demultiplexing reads results------"
while [ ! -e ${PID}_Demultiplexing.ok ]; do sleep 60 ; done
echo "------ /Demultiplexing reads ------"
