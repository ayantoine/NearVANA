#! /bin/bash

datetime1=$(date +%s)

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

ARG=$1
source $ARG
source $CONF
source $DATA
SDIR=${GITDIR}/Workflow

USE_PAIREND="$(boolean "${PAIREND}")"
USE_METADATA="$(boolean "${METADATA}")"
USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"
USE_KEEPUNASSIGNED="$(boolean "${UNASSIGNED}")"

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
	VAR_R1_FILE="${VARNAME}[$ID_R1]"
	TARGET_FILES=${!VAR_R1_FILE}
	if [ "$USE_PAIREND" = true ] ; then
		VAR_R2_FILE="${VARNAME}[$ID_R2]"
		TARGET_FILES+=" and "${!VAR_R2_FILE}
	fi
	VAR_DODE="${VARNAME}[$ID_DODE]"
	
	echo "------ Create Kmer library ------"
	echo "python ${SDIR}/CreateKmerList.py -m ${!VAR_DODE} -o ${!VAR_DODE}.kmer.tsv -v ${VARNAME}"
	python ${SDIR}/CreateKmerList.py -m ${!VAR_DODE} -o ${!VAR_DODE}.kmer.tsv -v ${VARNAME}
	echo "------ /Create Kmer library ------"
	
	echo "------ Count sequences ------"
	NB_SEQ=$(($(sed -n '$=' ${PID}_${VARNAME}_R1.fastq)/4))
	echo ${NB_SEQ} " present in ${PID}_${VARNAME}_R1.fastq"
	echo "------ /Count sequences ------"
	
	mkdir ${PID}_${VARNAME}_Demultiplexing
	
	echo "------ Make assignation ------"
	echo "python ${SDIR}/QsubAssignation.py -v ${VARNAME} -a ${ARG} -s ${SDIR} -k ${!VAR_DODE}.kmer.tsv -d ${PID}_${VARNAME}_Demultiplexing -o QsubAssignation.sh -c ${CONF} -q ${NB_SEQ} -p ${PID}"
	python ${SDIR}/QsubAssignation.py -v ${VARNAME} -a ${ARG} -s ${SDIR} -k ${!VAR_DODE}.kmer.tsv -d ${PID}_${VARNAME}_Demultiplexing -o QsubAssignation.sh -c ${CONF} -q ${NB_SEQ} -p ${PID}
	cat ./QsubAssignation.sh
	bash ./QsubAssignation.sh ${ARG}
	echo "------ /Make assignation ------"
	
	sleep 60
	
	echo "------ Merge subdata ------"
	cat ${PID}_${VARNAME}_Demultiplexing/*_Hyper_Identified* > ${PID}_${VARNAME}_Hyper_Identified.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Hypo_1_Identified* > ${PID}_${VARNAME}_Hypo_1_Identified.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Hypo_2_Identified* > ${PID}_${VARNAME}_Hypo_2_Identified.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Ambiguous_1* > ${PID}_${VARNAME}_Ambiguous_1.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Ambiguous_2* > ${PID}_${VARNAME}_Ambiguous_2.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Unidentified* > ${PID}_${VARNAME}_Unidentified.tsv
	echo "------ /Merge subdata ------"
	
	rm -r ${PID}_${VARNAME}_Demultiplexing/
	rm ./QsubAssignation.sh
	
	echo "------ Write output ------"
	echo "python ${SDIR}/ConcatenateFile.py -o ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -l ${PID}_${VARNAME}_Hyper_Identified.tsv,${PID}_${VARNAME}_Hypo_2_Identified.tsv,${PID}_${VARNAME}_Ambiguous_2.tsv,${PID}_${VARNAME}_Unidentified.tsv"
	python ${SDIR}/ConcatenateFile.py -o ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -l ${PID}_${VARNAME}_Hyper_Identified.tsv,${PID}_${VARNAME}_Hypo_2_Identified.tsv,${PID}_${VARNAME}_Ambiguous_2.tsv,${PID}_${VARNAME}_Unidentified.tsv
	echo "python ${SDIR}/ConcatenateFile.py -o ${PID}_${VARNAME}_Demultiplexing_Global.tsv -l ${PID}_${VARNAME}_Hyper_Identified.tsv,${PID}_${VARNAME}_Hypo_1_Identified.tsv,${PID}_${VARNAME}_Ambiguous_1.tsv,${PID}_${VARNAME}_Unidentified.tsv"
	python ${SDIR}/ConcatenateFile.py -o ${PID}_${VARNAME}_Demultiplexing_Global.tsv -l ${PID}_${VARNAME}_Hyper_Identified.tsv,${PID}_${VARNAME}_Hypo_1_Identified.tsv,${PID}_${VARNAME}_Ambiguous_1.tsv,${PID}_${VARNAME}_Unidentified.tsv
	echo "python ${SDIR}/CountDistribution.py -i ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -o ${PID}_${VARNAME}_Demultiplexing_Hyper_Distribution.tsv"
	python ${SDIR}/CountDistribution.py -i ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -o ${PID}_${VARNAME}_Demultiplexing_Hyper_Distribution.tsv
	echo "python ${SDIR}/CountDistribution.py -i ${PID}_${VARNAME}_Demultiplexing_Global.tsv -o ${PID}_${VARNAME}_Demultiplexing_Global_Distribution.tsv"
	python ${SDIR}/CountDistribution.py -i ${PID}_${VARNAME}_Demultiplexing_Global.tsv -o ${PID}_${VARNAME}_Demultiplexing_Global_Distribution.tsv
	echo "------ /Write output ------"
	
	echo "------ Bilan ------"
	wc -l ${PID}_${VARNAME}_Hyper_Identified.tsv ${PID}_${VARNAME}_Hypo_2_Identified.tsv ${PID}_${VARNAME}_Ambiguous_2.tsv ${PID}_${VARNAME}_Unidentified.tsv ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv | head -n -1
	echo $(expr $(cat ${PID}_${VARNAME}_R1.fastq ${PID}_${VARNAME}_R1.fastq | wc -l | cut -d " " -f1) / 4 )" sequences in "${TARGET_FILES}
	echo "------ /Bilan ------"

	echo "------ Store supplementary data ------"
	gzip -f ${PID}_${VARNAME}_Hypo_1_Identified.tsv > ${PID}_${VARNAME}_Hypo_1_Identified.tsv.gz
	gzip -f ${PID}_${VARNAME}_Hypo_2_Identified.tsv > ${PID}_${VARNAME}_Hypo_2_Identified.tsv.gz
	gzip -f ${PID}_${VARNAME}_Ambiguous_1.tsv > ${PID}_${VARNAME}_Ambiguous_1.tsv.gz
	gzip -f ${PID}_${VARNAME}_Ambiguous_2.tsv > ${PID}_${VARNAME}_Ambiguous_2.tsv.gz
	gzip -f ${PID}_${VARNAME}_Unidentified.tsv > ${PID}_${VARNAME}_Unidentified.tsv.gz
	gzip -f ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv > ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv.gz
	echo "------ /Store supplementary data ------"
	
	echo "------ Store unused data ------"
	if [ "$USE_KEEPUNASSIGNED" = true ] ; then
		gzip -f ${PID}_${VARNAME}_Hyper_Identified.tsv > ${PID}_${VARNAME}_Hyper_Identified.tsv.gz
	else
		gzip -f ${PID}_${VARNAME}_Demultiplexing_Global.tsv > ${PID}_${VARNAME}_Demultiplexing_Global.tsv.gz
	fi
	echo "------ /Store unused data ------"
	
done

> ${PID}.Demultiplexing.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Demultiplexing: "$delta > Time02.txt

