#! /bin/bash

set -e

err_report() {
    echo "Error on line $1"
}

trap 'err_report $LINENO' ERR

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

nb_jobs=$(ls ${PID}_BlastX/*_2.tab | wc -l)

echo "DEFINITION" > DEFINITION.txt
echo "ACCESSION" > ACCESSION.txt

if [ ! -d ${PID}"_log_Taxo" ] ; then mkdir ${PID}"_log_Taxo" ; fi
if [ ! -d "Taxo_Ok" ] ; then mkdir "Taxo_Ok" ; fi
touch ${PID}_nucleotide_TempDefDb.txt
touch ${PID}_protein_TempDefDb.txt
echo "$SCALL $SPARAM $SRENAME ${PID}_Taxo ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Taxo"/${PID}_Blast${TASK}.e${SPSEUDOTASKID} -o ${PID}"_log_Taxo"/${PID}_Blast${TASK}.o${SPSEUDOTASKID} ${SDIR}/GetTaxonomy.sh $ARG ${TASK}"
$SCALL $SPARAM $SRENAME ${PID}_Taxo ${STASKARRAY}1-${nb_jobs}${SMAXTASK}${SMAXSIMJOB} -e ${PID}"_log_Taxo"/${PID}_Blast${TASK}.e${SPSEUDOTASKID} -o ${PID}"_log_Taxo"/${PID}_Blast${TASK}.o${SPSEUDOTASKID} ${SDIR}/GetTaxonomy.sh $ARG ${TASK}
while true ; do
	if [ $(ls Taxo_Ok/ | wc -l) -eq 0 ]
		then
		nbr_ok=0
	else
		nbr_ok=$(ls Taxo_Ok/*_Taxo.ok | wc -l)
	fi
	if [ "${nbr_ok}" -eq "${nb_jobs}" ]
		then
		rm -r Taxo${TASK}_Ok
		break
	fi
	sleep 60
done

cat ${PID}_nucleotide_TempDefDb.txt >> ${NUCDEF}
cat ${PID}_protein_TempDefDb.txt >> ${PRODEF}

rm ${PID}_nucleotide_TempDefDb.txt
rm ${PID}_protein_TempDefDb.txt
rm DEFINITION.txt
rm ACCESSION.txt

touch ${PID}.Taxonomy.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Taxonomy: "$delta > Time10.txt
