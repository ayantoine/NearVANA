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

USE_PAIREND="$(boolean "${PAIREND}")"
USE_METADATA="$(boolean "${METADATA}")"
USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"
USE_KEEPUNASSIGNED="$(boolean "${UNASSIGNED}")"
USE_SUBSTRACTION="$(boolean "${SUBSTRACTION}")"
BLASTN="$(boolean "${BLASTN}")"
BLASTX="$(boolean "${BLASTX}")"
DIAMOND="$(boolean "${DIAMOND}")"
PREFILTER="$(boolean "${PREFILTER}")"

LIST_FILE=($DATA $ADAP $NUCACC $NUCDEF $PROACC $PRODEF $DBLINEAGE $VIRMINLEN $CONF)
if [ "$USE_SUBSTRACTION" = true ] ; then
	LIST_FILE+=($SUBS)
fi
for i in "${LIST_FILE[@]}"; do
	if [ ! -f $i ]; then
		echo "File $i does not exists"
		exit 1
	fi
done
echo "All referenced files in exists"
echo "------ check data"
source $DATA
NB_ITEM=1
ID_R1=0
if [ "$USE_PAIREND" = true ] ; then
	ID_R2=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi
if [ "$USE_MULTIPLEX" = true ] ; then
	ID_DODE=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
fi
if [ "$USE_METADATA" = true ] ; then
	ID_META=$NB_ITEM
	NB_ITEM=$((NB_ITEM+1))
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
		echo "With option MULTIPLEX set to $MULTIPLEX and PAIREND set to $PAIREND, $p must contains $NB_ITEM elements"
		exit 1
	fi
done
echo "All referenced files in $DATA exists"
echo "------ /Check Input existence ------"

echo "------ Show variable value ------"
echo "> Args details"
List_NONFILE=(PID DATA ADAP PAIREND METADATA MULTIPLEX UNASSIGNED SUBSTRACTION SUBS CONF BLASTX BLASTN DIAMOND PREFILTER)
for i in "${List_NONFILE[@]}"; do
	echo "$i: ${!i}"
done

echo
echo '> Conf details (/!\ beware bash interpretation for STASKID)'
source $CONF

SDIR=${GITDIR}/Analysis

echo "------ Check Input existence ------"
if [ ! -d $SDIR ] ; then
	echo "Directory SDIR ${SDIR} does not exists"
	exit 1
fi

LIST_PARAM=(NUCACC NUCDEF PROACC PRODEF DBLINEAGE VIRMINLEN SCALL SPARAM MULTICPU MULTIMEMORY SPARAM_MULTICPU STASKARRAY SMAXTASK SRENAME SMAXSIMJOB SMAXARRAYSIZE STASKID SPSEUDOTASKID VIRNTDB ALLNTDB  VIRPTDB ALLPTDB VIRPTDB_DIAMOND ALLPTDB_DIAMOND)
for i in "${LIST_PARAM[@]}"; do
	echo "$i: ${!i}"
done
echo "------ /Show variable value ------"

echo "------ Get Sample list ------"
declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	VAR_SAMPLE_FILE="${VARNAME}[$ID_DODE]"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done
if [ "$USE_KEEPUNASSIGNED" = true ] ; then
	SAMPLE_LIST+=("UnassignedReads")
fi
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
				echo "$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_${VARNAME}_R1_SplitReads.e -o Run_${VARNAME}_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_${VARNAME}_R1.fastq 1 ${VARNAME}"
				$SCALL $SPARAM $SRENAME ${PID}_R1_Run_SplitReads -e Run_${VARNAME}_R1_SplitReads.e -o Run_${VARNAME}_R1_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_${VARNAME}_R1.fastq 1 ${VARNAME}
				while [ ! -e ${PID}_${VARNAME}_R1.fastq.split.ok ]; do sleep 60 ; done
				rm ${PID}_${VARNAME}_R1.fastq
				if [ "$USE_PAIREND" = true ] ; then
					echo "$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_${VARNAME}_R2_SplitReads.e -o Run_${VARNAME}_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_${VARNAME}_R2.fastq 2 ${VARNAME}"
					$SCALL $SPARAM $SRENAME ${PID}_R2_Run_SplitReads -e Run_${VARNAME}_R2_SplitReads.e -o Run_${VARNAME}_R2_SplitReads.o ${SDIR}/Run_SplitReads.sh $ARG ${PID}_${VARNAME}_R2.fastq 2 ${VARNAME}
					while [ ! -e ${PID}_${VARNAME}_R2.fastq.split.ok ]; do sleep 60 ; done
					rm ${PID}_${VARNAME}_R2.fastq
				fi
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
	
	if [ "$USE_SUBSTRACTION" = true ] ; then
		if [ ! -f ${PID}.Substraction.ok ]; then
			echo -e "\t- PhiX Substraction : Susbract "${SUBS}
			echo "$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Substraction -e Substraction.e -o Substraction.o ${SDIR}/Substraction.sh $ARG"
			$SCALL $SPARAM_MULTICPU $SRENAME ${PID}_Substraction -e Substraction.e -o Substraction.o ${SDIR}/Substraction.sh $ARG
			while [ ! -e ${PID}.Substraction.ok ]; do sleep 60 ; done
			rm ${PID}_R1.Unsubstracted.fastq ${PID}_R2.Unsubstracted.fastq ${PID}_R0.Unsubstracted.fastq
			touch ${PID}.Substraction.ok
		else
			echo -e "\t- ${PID}.Substraction.ok existing, pas"
		fi	
	else
		mv ${PID}_R1.Unsubstracted.fastq ${PID}_R1.Substracted.fastq
		mv ${PID}_R2.Unsubstracted.fastq ${PID}_R2.Substracted.fastq
		mv ${PID}_R0.Unsubstracted.fastq ${PID}_R0.Substracted.fastq
	fi
	
else
	echo "${PID}.Cleaning.ok already existing, pass"
fi
echo "------ /Cleaning reads ------"

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

echo "------ Launch Blast/Diamond treatment------"
if [ "$BLASTN" = true ] ; then
	if [ ! -f ${PID}.BlastTreatmentN.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_BlastNTreatment -e Run_BlastNTreatment.e -o Run_BlastNTreatment.o ${SDIR}/Run_BlastTreatment.sh $ARG N"
		$SCALL $SPARAM $SRENAME ${PID}_BlastNTreatment -e Run_BlastNTreatment.e -o Run_BlastNTreatment.o ${SDIR}/Run_BlastTreatment.sh $ARG N
	else
		echo "${PID}.BlastTreatmentN.ok already existing, pass"
	fi
fi
if [ "$BLASTX" = true ] ; then
	if [ ! -f ${PID}.BlastTreatmentX.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_BlastXTreatment -e Run_BlastXTreatment.e -o Run_BlastXTreatment.o ${SDIR}/Run_BlastTreatment.sh $ARG X"
		$SCALL $SPARAM $SRENAME ${PID}_BlastXTreatment -e Run_BlastXTreatment.e -o Run_BlastXTreatment.o ${SDIR}/Run_BlastTreatment.sh $ARG X
	else
		echo "${PID}.BlastTreatmentX.ok already existing, pass"
	fi
fi
if [ "$DIAMOND" = true ] ; then
	if [ ! -f ${PID}.BlastTreatmentDiamond.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_DiamondTreatment -e Run_DiamondTreatment.e -o Run_DiamondTreatment.o ${SDIR}/Run_BlastTreatment.sh $ARG D"
		$SCALL $SPARAM $SRENAME ${PID}_DiamondTreatment -e Run_DiamondTreatment.e -o Run_DiamondTreatment.o ${SDIR}/Run_BlastTreatment.sh $ARG D
	else
		echo "${PID}.BlastTreatmentDiamond.ok already existing, pass"
	fi
fi
echo "------ /Launch Blast/Diamond treatment------"

if [ "$BLASTN" = true ] ; then
	while [ ! -e ${PID}.BlastTreatmentN.ok ]; do sleep 60 ; done
fi
if [ "$BLASTX" = true ] ; then
	while [ ! -e ${PID}.BlastTreatmentX.ok ]; do sleep 60 ; done
fi
if [ "$DIAMOND" = true ] ; then
	while [ ! -e ${PID}.BlastTreatmentD.ok ]; do sleep 60 ; done
fi

echo "------ Remove ToBlast data------"
rm -r ${PID}_ToBlast
echo "------ /Remove ToBlast data------"

echo "------ Produce basic stat------"
if [ ! -f ${PID}_All.Megahit_reverseAssembly.tsv ]; then
	echo "\t- Decompressing reverseAssembly archive"
	zcat ${PID}_All.Megahit_reverseAssembly.tsv.gz > ${PID}_All.Megahit_reverseAssembly.tsv
fi

if [ "$USE_MULTIPLEX" = true ] ; then
	if [ "$PREFILTER" = true ] ; then
		if [ "$BLASTD" = true ] ; then
			if [ ! -f ${PID}.StatBlastD.ok ]; then
				echo "python ${SDIR}/Extract4Stat.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA}"
				python ${SDIR}/Extract4Stat.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA}
			else
				echo "${PID}.StatBlastD.ok already existing, pass"
			fi
		fi
		if [ "$BLASTX" = true ] ; then
			if [ ! -f ${PID}.StatBlastX.ok ]; then
				echo "python ${SDIR}/Extract4Stat.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA}"
				python ${SDIR}/Extract4Stat.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA}
			else
				echo "${PID}.StatBlastX.ok already existing, pass"
			fi
		fi
		if [ "$BLASTN" = true ] ; then
			if [ ! -f ${PID}.StatBlastN.ok ]; then
				echo "python ${SDIR}/Extract4Stat.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA}"
				python ${SDIR}/Extract4Stat.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA}
			else
				echo "${PID}.StatBlastN.ok already existing, pass"
			fi
		fi
	else
		if [ "$BLASTD" = true ] ; then
			if [ ! -f ${PID}.StatBlastD.ok ]; then
				echo "python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv"
				python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv
			else
				echo "${PID}.StatBlastD.ok already existing, pass"
			fi
		fi
		if [ "$BLASTX" = true ] ; then
			if [ ! -f ${PID}.StatBlastX.ok ]; then
				echo "python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv"
				python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv
			else
				echo "${PID}.StatBlastX.ok already existing, pass"
			fi
		fi
		if [ "$BLASTN" = true ] ; then
			if [ ! -f ${PID}.StatBlastN.ok ]; then
				echo "python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv"
				python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv
			else
				echo "${PID}.StatBlastN.ok already existing, pass"
			fi
		fi
	fi
else
	if [ "$BLASTD" = true ] ; then
		if [ ! -f ${PID}.StatBlastD.ok ]; then
			echo "python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv"
			python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv
		else
			echo "${PID}.StatBlastD.ok already existing, pass"
		fi
	fi
	if [ "$BLASTX" = true ] ; then
		if [ ! -f ${PID}.StatBlastX.ok ]; then
			echo "python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv"
			python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv
		else
			echo "${PID}.StatBlastX.ok already existing, pass"
		fi
	fi
	if [ "$BLASTN" = true ] ; then
		if [ ! -f ${PID}.StatBlastN.ok ]; then
			echo "python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv"
			python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv
		else
			echo "${PID}.StatBlastN.ok already existing, pass"
		fi
	fi
fi

rm ${PID}_All.Megahit_reverseAssembly.tsv
echo "------ /Produce basic stat------"