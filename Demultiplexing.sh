#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
source $DATA

for VARNAME in "${PLATE[@]}"; do
	VAR_R1_FILE="${VARNAME}[0]"
	VAR_R2_FILE="${VARNAME}[1]"
	VAR_DODE="${VARNAME}[3]"
	
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
	bash ./QsubAssignation.sh
	echo "------ /Make assignation ------"
	
	sleep 60
	
	echo "------ Merge subdata ------"
	cat ${PID}_${VARNAME}_Demultiplexing/*_Hyper_Identified* > ${PID}_${VARNAME}_Hyper_Identified.tsv
	#cat ${PID}_${VARNAME}_Demultiplexing/*_Hypo_1_Identified* > ${PID}_${VARNAME}_Hypo_1_Identified.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Hypo_2_Identified* > ${PID}_${VARNAME}_Hypo_2_Identified.tsv
	#cat ${PID}_${VARNAME}_Demultiplexing/*_Ambiguous_1* > ${PID}_${VARNAME}_Ambiguous_1.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Ambiguous_2* > ${PID}_${VARNAME}_Ambiguous_2.tsv
	cat ${PID}_${VARNAME}_Demultiplexing/*_Unidentified* > ${PID}_${VARNAME}_Unidentified.tsv
	echo "------ /Merge subdata ------"
	
	rm -r ${PID}_${VARNAME}_Demultiplexing/
	rm ./QsubAssignation.sh
	
	echo "------ Write output ------"
	echo "python ${SDIR}/ConcatenateFile.py -o ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -l ${PID}_${VARNAME}_Hyper_Identified.tsv,${PID}_${VARNAME}_Hypo_2_Identified.tsv,${PID}_${VARNAME}_Ambiguous_2.tsv,${PID}_${VARNAME}_Unidentified.tsv"
	python ${SDIR}/ConcatenateFile.py -o ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -l ${PID}_${VARNAME}_Hyper_Identified.tsv,${PID}_${VARNAME}_Hypo_2_Identified.tsv,${PID}_${VARNAME}_Ambiguous_2.tsv,${PID}_${VARNAME}_Unidentified.tsv
	echo "python ${SDIR}/CountDistribution.py -i ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -o ${PID}_${VARNAME}_Demultiplexing_Hyper_Distribution.tsv"
	python ${SDIR}/CountDistribution.py -i ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv -o ${PID}_${VARNAME}_Demultiplexing_Hyper_Distribution.tsv
	echo "------ /Write output ------"
	
	echo "------ Bilan ------"
	wc -l ${PID}_${VARNAME}_Hyper_Identified.tsv ${PID}_${VARNAME}_Hypo_2_Identified.tsv ${PID}_${VARNAME}_Ambiguous_2.tsv ${PID}_${VARNAME}_Unidentified.tsv ${PID}_${VARNAME}_Demultiplexing_Hyper.tsv | head -n -1
	echo $(expr $(cat ${PID}_${VARNAME}_R1.fastq ${PID}_${VARNAME}_R1.fastq | wc -l | cut -d " " -f1) / 4 )" sequences in "${!VAR_R1_FILE}" and "${!VAR_R2_FILE}
	echo "------ /Bilan ------"

	echo "------ Store supplementary data ------"
	#gzip -f ${PID}_${VARNAME}_Hypo_1_Identified.tsv > ${PID}_Hypo_1_Identified.tsv.gz
	gzip -f ${PID}_${VARNAME}_Hypo_2_Identified.tsv > ${PID}_${VARNAME}_Hypo_2_Identified.tsv.gz
	#gzip -f $${PID}_${VARNAME}_Ambiguous_1.tsv > ${PID}_Ambiguous_1.tsv.gz
	gzip -f ${PID}_${VARNAME}_Ambiguous_2.tsv > ${PID}_${VARNAME}_Ambiguous_2.tsv.gz
	gzip -f ${PID}_${VARNAME}_Unidentified.tsv > ${PID}_${VARNAME}_Unidentified.tsv.gz
	echo "------ /Store supplementary data ------"
	
	#rm ${PID}_${VARNAME}_Hypo_2_Identified.tsv ${PID}_${VARNAME}_Ambiguous_2.tsv ${PID}_${VARNAME}_Unidentified.tsv
	
done

touch ${PID}.Demultiplexing.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Demultiplexing: "$delta > Time02.txt

