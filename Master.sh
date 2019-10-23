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

echo "------ Get Sample list ------"
declare -a SAMPLE_LIST
while read c1 leftovers; do
	SAMPLE_LIST+=($c1)
done < ${DODE}
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Sample list ------"

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
	touch ${PID}_Demultiplexing.ok
fi
while [ ! -e ${PID}_Demultiplexing.ok ]; do sleep 60 ; done
echo "------ /Demultiplexing reads -----"

echo "------ Cleaning ------"
if [ ! -f ${PID}_R1.Cleaned.fastq ] || [ ! -f ${PID}_R2.Cleaned.fastq ]; then
	for sampleId in "${SAMPLE_LIST[@]}"; do
		echo $sampleId
		mkdir $sampleId	
	done
	echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq"
	$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq 1 # Input: ${PID}_R1.fastq ${PID}_Hyper_Identified.tab
	echo "$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq"
	$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq 2 # Input: ${PID}_R2.fastq ${PID}_Hyper_Identified.tab
	while [ ! -e ${PID}_R1.fastq.split.ok ]; do sleep 60 ; done
	while [ ! -e ${PID}_R2.fastq.split.ok ]; do sleep 60 ; done
	echo "\t- Merge files"
	touch ${PID}_R1.Cleaned.fastq
	touch ${PID}_R2.Cleaned.fastq
	for sampleId in "${!SAMPLE_LIST[@]}"; do
		cat ${sampleId}/${sampleId}_${PID}_R1.fastq.split >> ${PID}_R1.Cleaned.fastq
		cat ${sampleId}/${sampleId}_${PID}_R2.fastq.split >> ${PID}_R2.Cleaned.fastq
	done
	touch ${PID}_Cleaning.okq
else
	echo "${PID}_R1.Cleaned.fastq and ${PID}_R2.Cleaned.fastq already existing, pass"
	touch ${PID}_Cleaning.ok
fi
echo "------ /Cleaning ------"


#echo "------ Trim adapters ------"
#if [ ! -f ${PID}_R1.Trimmed.fastq ] || [ ! -f ${PID}_R2.Trimmed.fastq ]; then
	#echo $SCALL $SPARAM $SRENAME ${PID}_Cleaning -e Cleaning.e -o Cleaning.o ${SDIR}/Cleaning.sh $ARG
	#$SCALL $SPARAM $SRENAME ${PID}_Cleaning -e Cleaning.e -o Cleaning.o ${SDIR}/Cleaning.sh $ARG # Input: ${PID}_R1.fastq ${PID}_R2.fastq $MID $PID Output: ${PID}_R1.Trimmed.fastq ${PID}_R2.Trimmed.fastq
#else
	#echo "${PID}_R1.Trimmed.fastq and ${PID}_R2.Trimmed.fastq  already existing, pass"
#fi
#while [ ! -e ${PID}_Cleaning.ok ]; do sleep 60 ; done
#echo "------ /Trim adapters ------"

