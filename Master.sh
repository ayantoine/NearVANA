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
if [ ! -f ${PID}_Demultiplexing.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing_Illumina_pe_V5.e -o Demultiplexing_Illumina_pe_V5.o ${SDIR}/Demultiplexing_Illumina_pe_V5.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing.e -o Demultiplexing.o ${SDIR}/Demultiplexing.sh $ARG # Input: ${PID}_R1.fastq ${PID}_R2.fastq $MID $PID $SDIR # Output: ${PID}_Demultiplexing.tab ${PID}_Demultiplexing_Distribution.tab ${PID}_Hyper_Identified.tab ${PID}_Hypo_1_Identified.tab ${PID}_Hypo_2_Identified.tab ${PID}_Ambiguous.tab ${PID}_Unidentified.tab
else
	echo "${PID}_Demultiplexing.ok existing, pass"
fi
while [ ! -e ${PID}_Demultiplexing.ok ]; do sleep 60 ; done
echo "------ /Demultiplexing reads -----"

#echo "------ Cleaning linkers ------"
#if [ ! -f ${PID}_Cleaning.ok ]; then
	#for sampleId in "${SAMPLE_LIST[@]}"; do
		#mkdir $sampleId	
	#done
	#echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq 1"
	#$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq 1 # Input: ${PID}_R1.fastq ${PID}_Hyper_Identified.tab
	#echo "$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq 2"
	#$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq 2 # Input: ${PID}_R2.fastq ${PID}_Hyper_Identified.tab
	#while [ ! -e ${PID}_R1.fastq.split.ok ]; do sleep 60 ; done
	#while [ ! -e ${PID}_R2.fastq.split.ok ]; do sleep 60 ; done
	#rm ${PID}_R1.fastq.split.ok ${PID}_R2.fastq.split.ok
	#echo -e "\t- Merge files"
	#touch ${PID}_R1.Cleaned.fastq
	#touch ${PID}_R2.Cleaned.fastq
	#for sampleId in "${SAMPLE_LIST[@]}"; do
		#cat ${sampleId}/${sampleId}_${PID}_R1.fastq.split >> ${PID}_R1.Cleaned.fastq
		#cat ${sampleId}/${sampleId}_${PID}_R2.fastq.split >> ${PID}_R2.Cleaned.fastq
	#done
	#touch ${PID}_Cleaning.ok
	#rm ${PID}_R1.fastq ${PID}_R2.fastq
#else
	#echo "${PID}_Cleaning.ok already existing, pass"
#fi
#echo "------ /Cleaning linkers ------"

#echo "------ Trim adapters ------"
#if [ ! -f ${PID}_Trimming.ok ]; then
	#echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_CutAdapt -e Run_R1_CutAdapt.e -o Run_R1_CutAdapt.o ${SDIR}/Run_Cutadapt.sh $ARG ${PID}_R1.fastq"
	#$SCALL $SPARAM $SRENAME ${PID}_Run_CutAdapt -e Run_CutAdapt.e -o Run_CutAdapt.o ${SDIR}/Run_Cutadapt.sh $ARG
	#while [ ! -e ${PID}.fastq.trim.ok ]; do sleep 60 ; done
	#rm ${PID}.fastq.trim.ok
	#echo -e "\t- Merge files"
	#touch ${PID}_R1.Trimmed.fastq
	#touch ${PID}_R2.Trimmed.fastq
	#for sampleId in "${SAMPLE_LIST[@]}"; do
		#cat ${sampleId}/${sampleId}_${PID}_R1.fastq.split.trim >> ${PID}_R1.Trimmed.fastq
		#cat ${sampleId}/${sampleId}_${PID}_R2.fastq.split.trim >> ${PID}_R2.Trimmed.fastq
	#done
	#touch ${PID}_Trimming.ok
	#rm ${PID}_R1.Cleaned.fastq ${PID}_R2.Cleaned.fastq
#else
	#echo "${PID}_Trimming.ok already existing, pass"
#fi
#echo "------ /Trim adapters ------"

#echo "------ PhiX Substraction ------"
#if [ ! -f ${PID}_Substraction.ok ]; then
	#echo -e "\t- Do deinterlacing"
	#echo "$SCALL $SPARAM $SRENAME ${PID}_Run_Correction -e Run_Correction.e -o Run_Correction.o ${SDIR}/Run_Correction.sh $ARG"
	#$SCALL $SPARAM $SRENAME ${PID}_Run_RetrievePair -e Run_RetrievePair.e -o Run_RetrievePair.o ${SDIR}/Run_RetrievePair.sh $ARG
	#while [ ! -e ${PID}.fastq.deinterlaced.ok ]; do sleep 60 ; done
	#rm ${PID}.fastq.deinterlaced.ok
	#echo -e "\t- Merge files"
	#touch ${PID}_R1.Unsubstracted.fastq
	#touch ${PID}_R2.Unsubstracted.fastq
	#touch ${PID}_R0.Unsubstracted.fastq
	#for sampleId in "${SAMPLE_LIST[@]}"; do
		#cat ${sampleId}/${sampleId}_${PID}_R1.fastq.split.trim.deinterlaced >> ${PID}_R1.Unsubstracted.fastq
		#cat ${sampleId}/${sampleId}_${PID}_R2.fastq.split.trim.deinterlaced >> ${PID}_R2.Unsubstracted.fastq
		#cat ${sampleId}/${sampleId}_${PID}_R0.fastq.split.trim.deinterlaced >> ${PID}_R0.Unsubstracted.fastq
	#done
	#echo -e "\t- Substract PhiX"
	#bash ${SDIR}/Substraction.sh $ARG
	#touch ${PID}_Substraction.ok
	#rm ${PID}_R1.Trimmed.fastq ${PID}_R2.Trimmed.fastq
	#rm ${PID}_R1.Unsubstracted.fastq ${PID}_R2.Unsubstracted.fastq ${PID}_R0.Unsubstracted.fastq
#else
	#echo "${PID}_Substraction.ok already existing, pass"
#fi
#echo "------ /PhiX Substraction------"

#echo "------ Reads correction ------"
#if [ ! -f ${PID}_Correction.ok ]; then
	#bash ${SDIR}/Correction.sh $ARG
	#touch ${PID}_Correction.ok
	#rm ${PID}_R1.Substracted.fastq ${PID}_R2.Substracted.fastq ${PID}_R0.Substracted.fastq
#else
	#echo "${PID}_Correction.ok already existing, pass"
#fi
#echo "------ /Reads correction------"













