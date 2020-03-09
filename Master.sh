#! /bin/bash

ARG=$1
#WARNING!! Source file is executed (Security, etc.)
source $ARG

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

BLASTN="$(boolean "${BLASTN}")"
BLASTX="$(boolean "${BLASTX}")"

if [ "$BLASTX" = false ] ; then
	if [ "$BLASTN" = false ] ; then
		echo "BlastX and BlastN set as FALSE, nothing to do\nexit"
		exit 1
	fi
fi

echo "------ Check Input existence ------"
if [ ! -d $SDIR ] ; then
	echo "Directory SDIR ${SDIR} does not exists"
	exit 1
fi
LIST_FILE=($R1 $R2 $ADAP $DODE $META $SUBS $NUCACC $NUCDEF $PROACC $PRODEF $DBLINEAGE $VIRMINLEN $CONF)
for i in "${LIST_FILE[@]}"; do
	if [ ! -f $i ]; then
		echo "File $i does not exists"
		exit 1
	fi
done

echo "> Args details"
List_NONFILE=(PID R1 R2 ADAP DODE META SUBS NUCACC NUCDEF PROACC PRODEF DBLINEAGE VIRMINLEN CONF)
for i in "${List_NONFILE[@]}"; do
	echo "$i: ${!i}"
done

echo
echo '> Conf details (/!\ beware bash interpretation for STASKID)'
source $CONF
LIST_PARAM=(SCALL SPARAM MULTICPU SPARAM_MULTICPU STASKARRAY SMAXTASK SRENAME SMAXSIMJOB SMAXARRAYSIZE STASKID SPSEUDOTASKID VIRNTDB ALLNTDB VIRPTDB ALLPTDB)
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
if [ ! -f ${PID}.extraction.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Extraction -e Extraction.e -o Extraction.o ${SDIR}/Gz_extraction.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Extraction -e Extraction.e -o Extraction.o ${SDIR}/Gz_extraction.sh $ARG # Output: ${PID}_R1.fastq ${PID}_R2.fastq
	while [ ! -e ${PID}.extraction.ok ]; do sleep 60 ; done
else
	echo "${PID}.extraction.ok already existing, pass"
fi
echo "------ /Extract .gz ------"

echo "------ Demultiplexing reads ------"
if [ ! -f ${PID}.Demultiplexing.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing_Illumina_pe_V5.e -o Demultiplexing_Illumina_pe_V5.o ${SDIR}/Demultiplexing_Illumina_pe_V5.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing.e -o Demultiplexing.o ${SDIR}/Demultiplexing.sh $ARG # Input: ${PID}_R1.fastq ${PID}_R2.fastq $MID $PID $SDIR # Output: ${PID}_Demultiplexing.tab ${PID}_Demultiplexing_Distribution.tab ${PID}_Hyper_Identified.tab ${PID}_Hypo_1_Identified.tab ${PID}_Hypo_2_Identified.tab ${PID}_Ambiguous.tab ${PID}_Unidentified.tab
	while [ ! -e ${PID}.Demultiplexing.ok ]; do sleep 60 ; done
else
	echo "${PID}.Demultiplexing.ok existing, pass"
fi
echo "------ /Demultiplexing reads -----"

echo "------ Cleaning reads ------"
if [ ! -f ${PID}.Cleaning.ok ]; then
	if [ ! -f ${PID}.SplitReads.ok ]; then
		echo -e "\t- Cleaning linkers"
		for sampleId in "${SAMPLE_LIST[@]}"; do
			mkdir $sampleId	
		done
		echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq 1"
		$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq 1 # Input: ${PID}_R1.fastq ${PID}_Hyper_Identified.tab
		echo "$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq 2"
		$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq 2 # Input: ${PID}_R2.fastq ${PID}_Hyper_Identified.tab
		while [ ! -e ${PID}_R1.fastq.split.ok ]; do sleep 60 ; done
		while [ ! -e ${PID}_R2.fastq.split.ok ]; do sleep 60 ; done
		rm ${PID}_R1.fastq ${PID}_R2.fastq
		touch ${PID}.SplitReads.ok
	else
		echo -e "\t- ${PID}.SplitReads.ok existing, pas"
	fi
	
	if [ ! -f ${PID}.CutAdapt.ok ]; then
		echo -e "\t- Trim adapters"
		echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_CutAdapt -e Run_R1_CutAdapt.e -o Run_R1_CutAdapt.o ${SDIR}/Run_Cutadapt.sh $ARG ${PID}_R1.fastq"
		$SCALL $SPARAM $SRENAME ${PID}_Run_CutAdapt -e Run_CutAdapt.e -o Run_CutAdapt.o ${SDIR}/Run_Cutadapt.sh $ARG
		while [ ! -e ${PID}.CutAdapt.ok ]; do sleep 60 ; done
		for sampleId in "${SAMPLE_LIST[@]}"; do
			rm ${sampleId}/${sampleId}_${PID}_R1.fastq.split
			rm ${sampleId}/${sampleId}_${PID}_R2.fastq.split
		done
	else
		echo -e "\t- ${PID}.CutAdapt.ok existing, pas"
	fi
	
	if [ ! -f ${PID}.Substraction-Deinterlacing.ok ]; then
		echo -e "\t- PhiX Substraction : deinterlacing"
		echo "$SCALL $SPARAM $SRENAME ${PID}_Run_Correction -e Run_RetrievePair.e -o Run_RetrievePair.o ${SDIR}/Run_RetrievePair.sh $ARG"
		$SCALL $SPARAM $SRENAME ${PID}_Run_RetrievePair -e Run_RetrievePair.e -o Run_RetrievePair.o ${SDIR}/Run_RetrievePair.sh $ARG
		while [ ! -e ${PID}.Deinterlacing.ok ]; do sleep 60 ; done
		rm ${PID}.Deinterlacing.ok
		
		for sampleId in "${SAMPLE_LIST[@]}"; do
			rm ${sampleId}/${sampleId}_${PID}_R1.fastq.split.trim
			rm ${sampleId}/${sampleId}_${PID}_R2.fastq.split.trim
		done
			
		echo -e "\t- PhiX Substraction : Merge deinterlaced subfiles"
		touch ${PID}_R1.Unsubstracted.fastq
		touch ${PID}_R2.Unsubstracted.fastq
		touch ${PID}_R0.Unsubstracted.fastq
		
		for sampleId in "${SAMPLE_LIST[@]}"; do
			cat ${sampleId}/${sampleId}_${PID}_R1.fastq.split.trim.deinterlaced >> ${PID}_R1.Unsubstracted.fastq
			cat ${sampleId}/${sampleId}_${PID}_R2.fastq.split.trim.deinterlaced >> ${PID}_R2.Unsubstracted.fastq
			cat ${sampleId}/${sampleId}_${PID}_R0.fastq.split.trim.deinterlaced >> ${PID}_R0.Unsubstracted.fastq
			rm -r ${sampleId}
		done
		touch ${PID}.Substraction-Deinterlacing.ok
	else
		echo -e "\t- ${PID}.Substraction-Deinterlacing.ok existing, pas"
	fi
	
	if [ ! -f ${PID}.Cleaning.ok ]; then
		echo -e "\t- PhiX Substraction : Susbract "${SUBS}
		echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Substraction -e Substraction.e -o Substraction.o ${SDIR}/Substraction.sh $ARG"
		$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Substraction -e Substraction.e -o Substraction.o ${SDIR}/Substraction.sh $ARG
		while [ ! -e ${PID}.Substraction.ok ]; do sleep 60 ; done
		rm ${PID}_R1.Unsubstracted.fastq ${PID}_R2.Unsubstracted.fastq ${PID}_R0.Unsubstracted.fastq
		touch ${PID}.Cleaning.ok
	else
		echo -e "\t- ${PID}.Cleaning.ok existing, pas"
	fi
else
	echo "${PID}.Cleaning.ok already existing, pass"
fi
echo "------ /Cleaning reads ------"

echo "------ Reads correction ------"
if [ ! -f ${PID}.Correction.ok ]; then
	echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Correction -e Correction.e -o Correction.o ${SDIR}/Correction.sh $ARG"
	$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Correction -e Correction.e -o Correction.o ${SDIR}/Correction.sh $ARG
	while [ ! -e ${PID}.Correction.ok ]; do sleep 60 ; done
	rm ${PID}_R1.Substracted.fastq ${PID}_R2.Substracted.fastq ${PID}_R0.Substracted.fastq
else
	echo "${PID}.Correction.ok already existing, pass"
fi
echo "------ /Reads correction------"

echo "------ Reads assembly ------"
if [ ! -f ${PID}.Assembly.ok ]; then
	echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Assembly -e Assembly.e -o Assembly.o ${SDIR}/Assembly.sh $ARG"
	$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Assembly -e Assembly.e -o Assembly.o ${SDIR}/Assembly.sh $ARG
	while [ ! -e ${PID}.Assembly.ok ]; do sleep 60 ; done
else
	echo "${PID}.Assembly.ok already existing, pass"
fi
echo "------ /Reads assembly------"

echo "------ Split fasta for Blast ------"
if [ ! -f ${PID}.SplitFasta.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_SplitFasta -e SplitFasta.e -o SplitFasta.o ${SDIR}/SplitFasta.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_SplitFasta -e SplitFasta.e -o SplitFasta.o ${SDIR}/SplitFasta.sh $ARG
	while [ ! -e ${PID}.SplitFasta.ok ]; do sleep 60 ; done
else
	echo "${PID}.SplitFasta.ok already existing, pass"
fi
echo "------ /Split fasta for Blast ------"

echo "------ Launch Blast treatment------"
if [ "$BLASTN" = true ] ; then
	if [ ! -f ${PID}.BlastN.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_BlastNTreatment -e Run_BlastNTreatment.e -o Run_BlastNTreatment.o ${SDIR}/Run_BlastNTreatment.sh $ARG N"
		$SCALL $SPARAM $SRENAME ${PID}_BlastNTreatment -e Run_BlastNTreatment.e -o Run_BlastNTreatment.o ${SDIR}/Run_BlastNTreatment.sh $ARG N
	else
		echo "${PID}.BlastN.ok already existing, pass"
	fi
fi
if [ "$BLASTX" = true ] ; then
	if [ ! -f ${PID}.BlastX.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_BlastXTreatment -e Run_BlastXTreatment.e -o Run_BlastXTreatment.o ${SDIR}/Run_BlastXTreatment.sh $ARG X"
		$SCALL $SPARAM $SRENAME ${PID}_BlastXTreatment -e Run_BlastXTreatment.e -o Run_BlastXTreatment.o ${SDIR}/Run_BlastXTreatment.sh $ARG X
	else
		echo "${PID}.BlastX.ok already existing, pass"
	fi
fi
rm -r ${PID}_ToBlast DEFINITION.txt ACCESSION.txt
echo "------ /Launch Blast treatment------"


#################################
#exit
#echo "$SCALL $SPARAM $SRENAME ${PID}_NTable -e Creation_TableAll.e -o Creation_TableAll.o ${SDIR}/CreateTableFusion.sh $ARG ${nb_jobs}"
#$SCALL $SPARAM $SRENAME ${PID}_AllTable -e Creation_TableAll.e -o Creation_TableAll.o ${SDIR}/CreateTableFusion.sh $ARG ${nb_jobs}
#perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_BlastAll_results.tab ${PID}_BlastAll_results.xlsx $((${#PID}+7))

##echo "------ Clean workdir ------"
##if [ ! -f ${PID}.Clean.ok ]; then
	##echo "$SCALL $SPARAM $SRENAME ${PID}_Clean -e CleanWorkDir.e -o CleanWorkDir.o ${SDIR}/CleanWorkDir.sh $ARG"
	## $SCALL $SPARAM $SRENAME ${PID}_Clean -e CleanWorkDir.e -o CleanWorkDir.o ${SDIR}/CleanWorkDir.sh $ARG
	##while [ ! -e ${PID}.Clean.ok ]; do sleep 60 ; done
##else
	##echo "${PID}.Clean.ok already existing, pass"
##fi
##echo "------ /Clean workdir ------"





