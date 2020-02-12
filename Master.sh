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

USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"

echo "------ Check Input existence ------"
if [ ! -d $SDIR ] ; then
	echo "Directory SDIR ${SDIR} does not exists"
	exit 1
fi
LIST_FILE=($DATA $ADAP $SUBS $NUCACC $NUCDEF $PROACC $PRODEF $DBLINEAGE $VIRMINLEN $CONF)
for i in "${LIST_FILE[@]}"; do
	if [ ! -f $i ]; then
		echo "File $i does not exists"
		exit 1
	fi
done
echo "All referenced files in $ARG exists"
echo "------ check data"
source $DATA
if [ "$USE_MULTIPLEX" = true ] ; then
	NB_ITEM=4
else
	NB_ITEM=2
fi

for VARNAME in "${PLATE[@]}"; do
	VAR=$VARNAME[@]
	#echo ${!VAR}
	eval "LEN=\${#$VAR[@]}"
	#echo "$len"
	if [ "$LEN" -eq "${NB_ITEM}" ]; then
		for i in "${p[@]}"; do
			if [ ! -f $i ]; then
				echo "File $i does not exists"
				exit 1
			fi
		done
	else
		echo "With option MULTIPLEX set to $MULTIPLEX, $p must contains $NB_ITEM elements"
		exit 1
	fi
done
echo "All referenced files in $DATA exists"
echo "------ /Check Input existence ------"

echo "------ Show variable value ------"
echo "> Args details"
List_NONFILE=(PID DATA MULTIPLEX ADAP SUBS NUCACC NUCDEF PROACC PRODEF DBLINEAGE VIRMINLEN CONF)
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
echo "------ /Show variable value ------"

echo "------ Get Sample list ------"
declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	VAR_SAMPLE_FILE="${VARNAME}[3]"
	#echo "${!VAR_SAMPLE_FILE}"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done
echo "${SAMPLE_LIST[@]}"
echo "------ /Get Sample list ------"


echo "------ Extract .gz ------"
if [ ! -f ${PID}.extraction.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Extraction -e Extraction.e -o Extraction.o ${SDIR}/Gz_extraction.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Extraction -e Extraction.e -o Extraction.o ${SDIR}/Gz_extraction.sh $ARG
	while [ ! -e ${PID}.extraction.ok ]; do sleep 60 ; done
else
	echo "${PID}.extraction.ok already existing, pass"
fi
echo "------ /Extract .gz ------"


if [ "$USE_MULTIPLEX" = true ] ; then
	echo "------ Demultiplexing reads ------"
	if [ ! -f ${PID}.Demultiplexing.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing.e -o Demultiplexing.o ${SDIR}/Demultiplexing.sh $ARG"
		$SCALL $SPARAM $SRENAME ${PID}_Demultiplexing -e Demultiplexing.e -o Demultiplexing.o ${SDIR}/Demultiplexing.sh $ARG
		while [ ! -e ${PID}.Demultiplexing.ok ]; do sleep 60 ; done
	else
		echo "${PID}.Demultiplexing.ok existing, pass"
	fi
	echo "------ /Demultiplexing reads -----"
fi

echo "------ Cleaning reads ------"
if [ ! -f ${PID}.Cleaning.ok ]; then
	if [ "$USE_MULTIPLEX" = true ] ; then
		if [ ! -f ${PID}.SplitReads.ok ]; then
			echo -e "\t- Cleaning linkers"
			for sampleId in "${SAMPLE_LIST[@]}"; do
				mkdir $sampleId	
			done
			for VARNAME in "${PLATE[@]}"; do
				echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_R1_SplitReads.e -o Run_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R1.fastq 1"
				$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_${VARNAME}_R1_SplitReads.e -o Run_${VARNAME}_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_${VARNAME}_R1.fastq 1 ${VARNAME}
				echo "$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_R2_SplitReads.e -o Run_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_R2.fastq 2"
				$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_${VARNAME}_R2_SplitReads.e -o Run_${VARNAME}_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_${VARNAME}_R2.fastq 2 ${VARNAME}
				while [ ! -e ${PID}_${VARNAME}_R1.fastq.split.ok ]; do sleep 60 ; done
				while [ ! -e ${PID}_${VARNAME}_R2.fastq.split.ok ]; do sleep 60 ; done
				rm ${PID}_${VARNAME}_R1.fastq ${PID}_${VARNAME}_R2.fastq
			done
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
	else
		if [ ! -f ${PID}.CutAdapt.ok ]; then
			echo -e "\t- Trim adapters"
			echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Run_CutAdapt_NM -e Run_CutAdapt.e -o Run_CutAdapt.o ${SDIR}/Run_Cutadapt_NM.sh $ARG"
			$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Run_CutAdapt_NM -e Run_CutAdapt.e -o Run_CutAdapt.o ${SDIR}/Run_Cutadapt_NM.sh $ARG
			while [ ! -e ${PID}.CutAdapt.ok ]; do sleep 60 ; done
		else
			echo -e "\t- ${PID}.CutAdapt.ok existing, pas"
		fi
		
		if [ ! -f ${PID}.Substraction-Deinterlacing.ok ]; then
			echo -e "\t- PhiX Substraction : deinterlacing"
			echo "$SCALL $SPARAM $SRENAME ${PID}_Run_RetrievePair_NM -e Run_RetrievePair.e -o Run_RetrievePair.o ${SDIR}/Run_RetrievePair_NM.sh $ARG"
			$SCALL $SPARAM $SRENAME ${PID}_Run_RetrievePair_NM -e Run_RetrievePair.e -o Run_RetrievePair.o ${SDIR}/Run_RetrievePair_NM.sh $ARG
			while [ ! -e ${PID}.Deinterlacing.ok ]; do sleep 60 ; done
			rm ${PID}.Deinterlacing.ok
			
			echo -e "\t- PhiX Substraction : Merge deinterlaced subfiles"
			touch ${PID}_R1.Unsubstracted.fastq
			touch ${PID}_R2.Unsubstracted.fastq
			touch ${PID}_R0.Unsubstracted.fastq
			
			for VARNAME in "${PLATE[@]}"; do
				cat ${PID}_${VARNAME}_R1.fastq.trim.deinterlaced >> ${PID}_R1.Unsubstracted.fastq
				cat ${PID}_${VARNAME}_R2.fastq.trim.deinterlaced >> ${PID}_R2.Unsubstracted.fastq
				cat ${PID}_${VARNAME}_R0.fastq.trim.deinterlaced >> ${PID}_R0.Unsubstracted.fastq
				rm ${PID}_${VARNAME}_R*.fastq.trim.deinterlaced
			done
			touch ${PID}.Substraction-Deinterlacing.ok
		else
			echo -e "\t- ${PID}.Substraction-Deinterlacing.ok existing, pas"
		fi	
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


########################################################################
exit

echo "------ Reads assembly ------"
if [ ! -f ${PID}.Assembly.ok ]; then
	echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Assembly -e Assembly.e -o Assembly.o ${SDIR}/Assembly.sh $ARG"
	$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Assembly -e Assembly.e -o Assembly.o ${SDIR}/Assembly.sh $ARG
	while [ ! -e ${PID}.Assembly.ok ]; do sleep 60 ; done
else
	echo "${PID}.Assembly.ok already existing, pass"
fi
echo "------ /Reads assembly------"

echo "------ Launch Blast ------"
if [ ! -f ${PID}.Blast.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Run_Blast -e Run_Blast.e -o Run_Blast.o ${SDIR}/Run_Blast.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Run_Blast -e Run_Blast.e -o Run_Blast.o ${SDIR}/Run_Blast.sh $ARG
	while [ ! -e ${PID}.Blast.ok ]; do sleep 60 ; done
else
	echo "${PID}.Blast.ok already existing, pass"
fi
echo "------ /Launch Blast------"

echo "------ Retrieve Taxonomy data ------"
if [ ! -f ${PID}.Taxonomy.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Run_Taxo -e Run_Taxo.e -o Run_Taxo.o ${SDIR}/Run_Taxo.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Run_Taxo -e Run_Taxo.e -o Run_Taxo.o ${SDIR}/Run_Taxo.sh $ARG
	while [ ! -e ${PID}.Taxonomy.ok ]; do sleep 60 ; done
else
	echo "${PID}.Taxonomy.ok already existing, pass"
fi
echo "------ /Retrieve Taxonomy data ------"

echo "------ Create table ------"
if [ ! -f ${PID}.Table.ok ]; then
	echo "$SCALL $SPARAM $SRENAME ${PID}_Table -e Table.e -o Table.o ${SDIR}/Run_CreateTable.sh $ARG"
	$SCALL $SPARAM $SRENAME ${PID}_Table -e Run_Table.e -o Run_Table.o ${SDIR}/Run_CreateTable.sh $ARG
	while [ ! -e ${PID}.Table.ok ]; do sleep 60 ; done
else
	echo "${PID}.Table.ok already existing, pass"
fi
echo "------ /Create table ------"

##echo "------ Clean workdir ------"
##if [ ! -f ${PID}.Clean.ok ]; then
	##echo "$SCALL $SPARAM $SRENAME ${PID}_Clean -e CleanWorkDir.e -o CleanWorkDir.o ${SDIR}/CleanWorkDir.sh $ARG"
	## $SCALL $SPARAM $SRENAME ${PID}_Clean -e CleanWorkDir.e -o CleanWorkDir.o ${SDIR}/CleanWorkDir.sh $ARG
	##while [ ! -e ${PID}.Clean.ok ]; do sleep 60 ; done
##else
	##echo "${PID}.Clean.ok already existing, pass"
##fi
##echo "------ /Clean workdir ------"





