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
SDIR=${GITDIR}/Workflow

source $DATA

USE_PAIREND="$(boolean "${PAIREND}")"
USE_METADATA="$(boolean "${METADATA}")"
USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"

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

USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"

echo "------ Megahit ------"
#time megahit --k-list 21,33,55,77,99 -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -r ${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${PID}"_log_Assembly-Megahit"
python ${SDIR}/Hack_SplitSubstractedBySample.py -i ${PID}_R1.Substracted.fastq
python ${SDIR}/Hack_SplitSubstractedBySample.py -i ${PID}_R2.Substracted.fastq
python ${SDIR}/Hack_SplitSubstractedBySample.py -i ${PID}_R0.Substracted.fastq

declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	echo "${VARNAME} ${ID_DODE}"
	VAR_SAMPLE_FILE="${VARNAME}[$ID_DODE]"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done

echo "DEBUG"
echo "${SAMPLE_LIST[@]}"


touch "Hack_"${PID}"_All.Megahit_rejectedContigs.fa"
touch "Hack_"${PID}"_All.Megahit_ambigousReads.tsv"
touch "Hack_"${PID}"_All.Megahit_unmappedReads.tsv"

touch "Hack_"${PID}"_All.Megahit_reverseAssembly.tsv"
touch "Hack_"${PID}"_All.Megahit_contigs.fa"
touch "Hack_"${PID}"_R1.Megahit_unassembled.fastq"
touch "Hack_"${PID}"_R2.Megahit_unassembled.fastq"
touch "Hack_"${PID}"_R0.Megahit_unassembled.fastq"


for sampleId in "${SAMPLE_LIST[@]}"; do
	echo "${sampleId}"
	
	echo "megahit --k-list 21,33,55,77,99 -1 ${sampleId}/${PID}_R1.Substracted.fastq -2 ${sampleId}/${PID}_R2.Substracted.fastq -r ${sampleId}/${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${sampleId}/${PID}_log_Assembly-Megahit"
	megahit --k-list 21,33,55,77,99 -1 ${sampleId}/${PID}_R1.Substracted.fastq -2 ${sampleId}/${PID}_R2.Substracted.fastq -r ${sampleId}/${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${sampleId}/${PID}"_log_Assembly-Megahit"
	
	#echo "python ${SDIR}/Hack_RenameFasta.py -i ${sampleId}/${PID}_log_Assembly-Megahit/final.contigs.fa -o ${sampleId}/${PID}_Temp.Megahit_contigs.fa -s ${sampleId}"
	#python ${SDIR}/Hack_RenameFasta.py -i ${sampleId}/${PID}_log_Assembly-Megahit/final.contigs.fa -o ${sampleId}/${PID}_Temp.Megahit_contigs.fa -s ${sampleId}
	
	echo "bowtie2-build --threads ${MULTICPU} ${sampleId}/${PID}_Temp.Megahit_contigs.fa ${sampleId}/${PID}_Temp.Megahit_contigs.fa"
	bowtie2-build --threads ${MULTICPU} ${sampleId}/${PID}_Temp.Megahit_contigs.fa ${sampleId}/${PID}_Temp.Megahit_contigs.fa
	if [ -s "${PID}_R0.Substracted.fastq" ]; then
		echo "bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${sampleId}/${PID}_Temp.Megahit_contigs.fa -1 ${sampleId}/${PID}_R1.Substracted.fastq -2 ${sampleId}/${PID}_R2.Substracted.fastq -U ${sampleId}/${PID}_R0.Substracted.fastq -S ${sampleId}/reads2contigs.sam"
		bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${sampleId}/${PID}_Temp.Megahit_contigs.fa -1 ${sampleId}/${PID}_R1.Substracted.fastq -2 ${sampleId}/${PID}_R2.Substracted.fastq -U ${sampleId}/${PID}_R0.Substracted.fastq -S ${sampleId}/reads2contigs.sam
	else
		echo "bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${sampleId}/${PID}_Temp.Megahit_contigs.fa -1 ${sampleId}/${PID}_R1.Substracted.fastq -2 ${sampleId}/${PID}_R2.Substracted.fastq -S ${sampleId}/reads2contigs.sam"
		bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${sampleId}/${PID}_Temp.Megahit_contigs.fa -1 ${sampleId}/${PID}_R1.Substracted.fastq -2 ${sampleId}/${PID}_R2.Substracted.fastq -S ${sampleId}/reads2contigs.sam
	fi
	rm ${sampleId}/*.bt2
	
	scp ${sampleId}/reads2contigs.sam ./reads2contigs.sam
	scp ${sampleId}/${PID}_Temp.Megahit_contigs.fa ./${PID}_Temp.Megahit_contigs.fa
	scp ${sampleId}/*.fastq ./
	
	echo "python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX} -s ${SampleID}"
	python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX} -s ${SampleID}
	
	cat ${PID}"_All.Megahit_rejectedContigs.fa" >> "Hack_"${PID}"_All.Megahit_rejectedContigs.fa"
	cat ${PID}"_All.Megahit_ambigousReads.tsv" >> "Hack_"${PID}"_All.Megahit_ambigousReads.tsv"
	cat ${PID}"_All.Megahit_unmappedReads.tsv" >> "Hack_"${PID}"_All.Megahit_unmappedReads.tsv"

	cat ${PID}"_All.Megahit_reverseAssembly.tsv" >> "Hack_"${PID}"_All.Megahit_reverseAssembly.tsv"
	cat ${PID}"_All.Megahit_contigs.fa" >> "Hack_"${PID}"_All.Megahit_contigs.fa"
	cat ${PID}"_R1.Megahit_unassembled.fastq" >> "Hack_"${PID}"_R1.Megahit_unassembled.fastq"
	cat ${PID}"_R2.Megahit_unassembled.fastq" >> "Hack_"${PID}"_R2.Megahit_unassembled.fastq"
	cat ${PID}"_R0.Megahit_unassembled.fastq" >> "Hack_"${PID}"_R0.Megahit_unassembled.fastq"
	
	echo "remove ${sampleId}"
	#rm -r ${sampleId}/
	
	rm ${PID}"_All.Megahit_rejectedContigs.fa" 
	rm ${PID}"_All.Megahit_ambigousReads.tsv"
	rm ${PID}"_All.Megahit_unmappedReads.tsv"
	
	rm ${PID}"_All.Megahit_reverseAssembly.tsv"
	rm ${PID}"_All.Megahit_contigs.fa"
	rm ${PID}"_R1.Megahit_unassembled.fastq"
	rm ${PID}"_R2.Megahit_unassembled.fastq"
	rm ${PID}"_R0.Megahit_unassembled.fastq"
	
done

rm ./reads2contigs.sam
rm ./${PID}_Temp.Megahit_contigs.fa

mv "Hack_"${PID}"_All.Megahit_rejectedContigs.fa" ${PID}"_All.Megahit_rejectedContigs.fa" 
mv "Hack_"${PID}"_All.Megahit_ambigousReads.tsv" ${PID}"_All.Megahit_ambigousReads.tsv"
mv "Hack_"${PID}"_All.Megahit_unmappedReads.tsv" ${PID}"_All.Megahit_unmappedReads.tsv"

mv "Hack_"${PID}"_All.Megahit_reverseAssembly.tsv" ${PID}"_All.Megahit_reverseAssembly.tsv"
mv "Hack_"${PID}"_All.Megahit_contigs.fa" ${PID}"_All.Megahit_contigs.fa"
mv "Hack_"${PID}"_R1.Megahit_unassembled.fastq" ${PID}"_R1.Megahit_unassembled.fastq"
mv "Hack_"${PID}"_R2.Megahit_unassembled.fastq" ${PID}"_R2.Megahit_unassembled.fastq"
mv "Hack_"${PID}"_R0.Megahit_unassembled.fastq" ${PID}"_R0.Megahit_unassembled.fastq"


echo "------ /Megahit ------"

if [ "$USE_MULTIPLEX" = true ] ; then
	if [ ! -f ${PID}.Stat_Assembly.ok ]; then
		echo "------ Write stat ------"
		echo "$SCALL $SPARAM_HEAVY $SRENAME ${PID}_${TASK}Stat -e Stat_Assembly.e -o Stat_Assembly.o ${SDIR}/CountAssemblyStat.sh $ARG"
		$SCALL $SPARAM_HEAVY $SRENAME ${PID}_${TASK}Stat -e Stat_Assembly.e -o Stat_Assembly.o ${SDIR}/CountAssemblyStat.sh $ARG
		while [ ! -e ${PID}.Stat_Assembly.ok ]; do sleep 60 ; done
		echo "------ /Write stat ------"
	fi
fi

echo "------ Compress Corrected.fastq ------"
gzip -f ${PID}_R0.Substracted.fastq
gzip -f ${PID}_R1.Substracted.fastq
gzip -f ${PID}_R2.Substracted.fastq
echo "------ /Compress Corrected.fastq ------"

echo "------ Merge Assembly ------"
mv ${PID}_All.Megahit_contigs.fa ${PID}_All.fa
echo "------ /Merge Assembly ------"

echo "------ Compress Megahit output ------"
gzip -f ${PID}_All.Megahit_reverseAssembly.tsv
gzip -f ${PID}_All.Megahit_ambigousReads.tsv
gzip -f ${PID}_All.Megahit_rejectedContigs.fa
gzip -f ${PID}_All.Megahit_unmappedReads.tsv
echo "------ /Compress Megahit output ------"

touch ${PID}.Assembly.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Assembly: "$delta > Time07.txt
