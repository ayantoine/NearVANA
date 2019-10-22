#! /bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

echo "------ Split fastq ------"
split -a 3 -d -l 4000000 ${PID}_R1.fastq 1Clean_ #1,000,000 sequences by subfiles
split -a 3 -d -l 4000000 ${PID}_R2.fastq 2Clean_
echo "------ /Split fastq ------"

mkdir ${PID}_Cleaning
mv 1Clean* ${PID}_Cleaning/
mv 2Clean* ${PID}_Cleaning/

NumberFile=$(ls ${PID}_Cleaning/ | wc -l)
NumberPairFile=$(($NumberFile/2))

echo "------ Parallelize CutAdapt ------"
if [ ! -d "CutAdapt_Ok" ] ; then mkdir "CutAdapt_Ok" ; fi
echo $SCALL $SPARAM $SRENAME ${PID}_Run_Cutadapt ${STASKARRAY}1-${NumberPairFile}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_Cleaning"/${PID}_Run_Cutadapt.e${SPSEUDOTASKID} -o ${PID}"_Cleaning"/${PID}_Run_Cutadapt.o${SPSEUDOTASKID} ${SDIR}/Run_Cutadapt.sh $ARG
$SCALL $SPARAM $SRENAME ${PID}_Run_Cutadapt ${STASKARRAY}1-${NumberPairFile}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_Cleaning"/${PID}_Run_Cutadapt.e${SPSEUDOTASKID} -o ${PID}"_Cleaning"/${PID}_Run_Cutadapt.o${SPSEUDOTASKID} ${SDIR}/Run_Cutadapt.sh $ARG
while true ; do
	if [ $(ls CutAdapt_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls CutAdapt_Ok/*_CutAdapat.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq "${NumberPairFile}" ]
		then
		rm -r CutAdapt_Ok
		break
	fi
	sleep 60
done
echo "------ /Parallelize CutAdapt ------"

echo "------ Merge subdata ------"
cat ${PID}_Cleaning/*R1.Trim1.fastq > ${PID}_R1.Trimmed.fastq
cat ${PID}_Cleaning/*R2.Trim1.fastq > ${PID}_R2.Trimmed.fastq
echo "------ /Merge subdata ------"

rm ${MID}.C.tab
rm ${PID}_TempC1.tab

rm ${PID}_Cleaning/*Clean_*
rm ${PID}_Cleaning/*Trim1*

touch ${PID}_Cleaning.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Cleaning: "$delta > Time03.txt
