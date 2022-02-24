#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"

echo "------ Megahit ------"
#time megahit --k-list 21,33,55,77,99 -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -r ${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${PID}"_log_Assembly-Megahit"
python ${SDIR}/Hack_SplitSubstractedBySample.py -i ${PID}_R1.Substracted.fastq
python ${SDIR}/Hack_SplitSubstractedBySample.py -i ${PID}_R2.Substracted.fastq
python ${SDIR}/Hack_SplitSubstractedBySample.py -i ${PID}_R0.Substracted.fastq

declare -a SAMPLE_LIST
for VARNAME in "${PLATE[@]}"; do
	VAR_SAMPLE_FILE="${VARNAME}[$ID_DODE]"
	while read c1 leftovers; do
		SAMPLE_LIST+=(${VARNAME}${c1})
	done < ${!VAR_SAMPLE_FILE}
done


touch ${PID}"_All.Megahit_rejectedContigs.fa"
touch ${PID}"_All.Megahit_ambigousReads.tsv"
touch ${PID}"_All.Megahit_unmappedReads.tsv"

touch ${PID}"_All.Megahit_reverseAssembly.tsv"
touch ${PID}"_All.Megahit_contigs.fa"
touch ${PID}"_R1.Megahit_unassembled.fastq"
touch ${PID}"_R2.Megahit_unassembled.fastq"
touch ${PID}"_R0.Megahit_unassembled.fastq"


for sampleId in "${SAMPLE_LIST[@]}"; do
	cd $sampleId/
	
	megahit --k-list 21,33,55,77,99 -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -r ${PID}_R0.Substracted.fastq -m ${MULTIMEMORY} -t ${MULTICPU} -o ${PID}"_log_Assembly-Megahit"
	
	#mv ${PID}_log_Assembly-Megahit/final.contigs.fa ${PID}_Temp.Megahit_contigs.fa
	python ${SDIR}/Hack_RenameFastq.py -i ${PID}_log_Assembly-Megahit/final.contigs.fa -o ${PID}_Temp.Megahit_contigs.fa -s $sampleId
	
	
	bowtie2-build --threads ${MULTICPU} ${PID}_Temp.Megahit_contigs.fa ${PID}_Temp.Megahit_contigs.fa
	if [ -s "${PID}_R0.Substracted.fastq" ]; then
		bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -U ${PID}_R0.Substracted.fastq -S reads2contigs.sam
	else
		bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -S reads2contigs.sam
	fi
	rm *.bt2
	echo "python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}"
	python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}
	
	cd ..
	
	cat $sampleId/${PID}"_All.Megahit_rejectedContigs.fa" >> ${PID}"_All.Megahit_rejectedContigs.fa"
	cat $sampleId/${PID}"_All.Megahit_ambigousReads.tsv" >> ${PID}"_All.Megahit_ambigousReads.tsv"
	cat $sampleId/${PID}"_All.Megahit_unmappedReads.tsv" >> ${PID}"_All.Megahit_unmappedReads.tsv"

	cat $sampleId/${PID}"_All.Megahit_reverseAssembly.tsv" >> ${PID}"_All.Megahit_reverseAssembly.tsv"
	cat $sampleId/${PID}"_All.Megahit_contigs.fa" >> ${PID}"_All.Megahit_contigs.fa"
	cat $sampleId/${PID}"_R1.Megahit_unassembled.fastq" >> ${PID}"_R1.Megahit_unassembled.fastq"
	cat $sampleId/${PID}"_R2.Megahit_unassembled.fastq" >> ${PID}"_R2.Megahit_unassembled.fastq"
	cat $sampleId/${PID}"_R0.Megahit_unassembled.fastq" >> ${PID}"_R0.Megahit_unassembled.fastq"

	rm -r $sampleId/

done


echo "------ /Megahit ------"

#mv ${PID}_log_Assembly-Megahit/final.contigs.fa ${PID}_Temp.Megahit_contigs.fa

#echo "------ Megahit reverse-mapping ------"
#bowtie2-build --threads ${MULTICPU} ${PID}_Temp.Megahit_contigs.fa ${PID}_Temp.Megahit_contigs.fa
#if [ -s "${PID}_R0.Substracted.fastq" ]; then
	#bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -U ${PID}_R0.Substracted.fastq -S reads2contigs.sam
#else
	#bowtie2 --threads ${MULTICPU} --end-to-end --very-sensitive -x ${PID}_Temp.Megahit_contigs.fa -1 ${PID}_R1.Substracted.fastq -2 ${PID}_R2.Substracted.fastq -S reads2contigs.sam
#fi
#rm *.bt2
#echo "python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}"
#python ${SDIR}/MappingReverseMegahit.py -p ${PID} -i reads2contigs.sam -m ${MULTIPLEX}
#echo "------ /Megahit reverse-mapping ------"

#rm ${PID}_Temp.Megahit_contigs.fa reads2contigs.sam

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
