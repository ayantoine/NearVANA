#! /bin/bash

ARG=$1
source $ARG
source $CONF

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}


USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"

#N or X
TASK=$2


if [ ${Task} == D ]; then
	echo "------ Launch Diamond by task ------"
	if [ ! -f {PID}.Blast${TASK}.ok ]; then
		echo ""
		$SCALL $SPARAM $SRENAME ${PID}_DiamondTask -e Run_DiamondTask.e -o Run_DiamondTask.o ${SDIR}/Run_DiamondTask.sh $ARG
		while [ ! -e ${PID}.Blast${TASK}.ok ]; do sleep 60 ; done
	else
		echo "{PID}.Blast${TASK}.ok already existing, pass"
	fi
	echo "------ /Launch Diamond by task ------"
else
	echo "------ Launch Blast by task ------"
	if [ ! -f {PID}.Blast${TASK}.ok ]; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_${TASK}_Blast -e Run_Blast${TASK}.e -o Run_Blast${TASK}.o ${SDIR}/Run_BlastTask.sh $ARG ${TASK}"
		$SCALL $SPARAM $SRENAME ${PID}_${TASK}_Blast -e Run_Blast${TASK}.e -o Run_Blast${TASK}.o ${SDIR}/Run_BlastTask.sh $ARG ${TASK}
		while [ ! -e ${PID}.Blast${TASK}.ok ]; do sleep 60 ; done
	else
		echo "{PID}.Blast${TASK}.ok already existing, pass"
	fi
	echo "------ /Launch Blast by task ------"
fi

echo "------ Retrieve Taxonomy data ------"
if [ ! -f ${PID}.Taxonomy${TASK}.ok ]; then
	nb_jobs=$(ls ${PID}_Blast${TASK}/*_2.tab | wc -l)

	echo "DEFINITION" > DEFINITION.txt
	echo "ACCESSION" > ACCESSION.txt

	if [ ! -d ${PID}"_log_Taxo${TASK}" ] ; then mkdir ${PID}"_log_Taxo${TASK}" ; fi
	if [ ! -d "Taxo${TASK}_Ok" ] ; then mkdir "Taxo${TASK}_Ok" ; fi
	touch ${PID}_${TASK}_TempDefDb.txt
	echo "$SCALL $SPARAM $SRENAME ${PID}_Taxo${TASK} ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Taxo${TASK}"/${PID}_Taxo${TASK}.e${SPSEUDOTASKID} -o ${PID}"_log_Taxo${TASK}"/${PID}_Taxo${TASK}.o${SPSEUDOTASKID} ${SDIR}/GetTaxonomy.sh $ARG ${TASK}"
	$SCALL $SPARAM $SRENAME ${PID}_Taxo${TASK} ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Taxo${TASK}"/${PID}_Taxo${TASK}.e${SPSEUDOTASKID} -o ${PID}"_log_Taxo${TASK}"/${PID}_Taxo${TASK}.o${SPSEUDOTASKID} ${SDIR}/GetTaxonomy.sh $ARG ${TASK}
	while true ; do
		if [ $(ls Taxo${TASK}_Ok/ | wc -l) -eq 0 ]
			then
			nbr_ok=0
		else
			nbr_ok=$(ls Taxo${TASK}_Ok/*_Taxo.ok | wc -l)
		fi
		if [ "${nbr_ok}" -eq "${nb_jobs}" ]
			then
			rm -r Taxo${TASK}_Ok
			break
		fi
		sleep 60
	done
	
	cat ${PID}_${TASK}_TempDefDb.txt >> ${NUCDEF}
	rm ${PID}_${TASK}_TempDefDb.txt
	
	touch ${PID}.Taxonomy${TASK}.ok
else
	echo "${PID}.Taxonomy${TASK}.ok already existing, pass"
fi
echo "------ /Retrieve Taxonomy data ------"

echo "------ Create table ------"
if [ ! -f ${PID}.creation${TASK}.ok ]; then
	if [ ! -f ${PID}_All.SPAdes.contigs2sample.tsv ]; then
		echo "------ Create short-list ------"
		cut -f2,4 ${PID}_All.SPAdes_reverseAssembly.tsv | sort -u > ${PID}_All.SPAdes.contigs2sample.tsv
		cut -f2,4 ${PID}_All.FLASH_reverseAssembly.tsv | sort -u > ${PID}_All.FLASH.contigs2sample.tsv
		echo "------ /Create short-list ------"
	fi
	datetime1=$(date +%s)
	if [ "$USE_MULTIPLEX" = true ] ; then
		echo "$SCALL $SPARAM $SRENAME ${PID}_${TASK}Table -e Creation_Table${TASK}.e -o Creation_Table${TASK}.o ${SDIR}/CreateTable.sh $ARG ${TASK} ${nb_jobs}"
		$SCALL $SPARAM $SRENAME ${PID}_${TASK}Table -e Creation_Table${TASK}.e -o Creation_Table${TASK}.o ${SDIR}/CreateTable.sh $ARG ${TASK} ${nb_jobs}
		while [ ! -e ${PID}.creation${TASK}.ok ]; do sleep 60 ; done
	else
		echo "$SCALL $SPARAM $SRENAME ${PID}_${TASK}Table -e Creation_Table${TASK}.e -o Creation_Table${TASK}.o ${SDIR}/CreateTable_NM.sh $ARG ${TASK} ${nb_jobs}"
		$SCALL $SPARAM $SRENAME ${PID}_${TASK}Table -e Creation_Table${TASK}.e -o Creation_Table${TASK}.o ${SDIR}/CreateTable_NM.sh $ARG ${TASK} ${nb_jobs}
		while [ ! -e ${PID}.creation${TASK}.ok ]; do sleep 60 ; done
	fi
	echo "------ Xlsx conversion ------"
	perl -I ${SDIR} ${SDIR}/Tab2Xls.pl ${PID}_Blast${TASK}_results.tab ${PID}_Blast${TASK}_results.xlsx $((${#PID}+7))
	echo "------ /Xlsx conversion ------"
	datetime2=$(date +%s)
	delta=$((datetime2 - datetime1))
	echo "Time Table: "$delta > Time12.txt
else
	echo "${PID}.creation${TASK}.ok already existing, pass"
fi
echo "------ /Create table ------"
