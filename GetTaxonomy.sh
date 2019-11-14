#!/bin/bash

ARG=$1
source $ARG
source $CONF

TASK_ARRAY=(N X)

for Task in ${TASK_ARRAY[@]}; do
    FilePath=${PID}_Blast${Task}/${PID}_All.fa.${STASKID}.Blast${Task}_2.tab
    if [ ${Task} == X ] ; then
	DBTARGET=${PROACC}
	DBDEF=${PRODEF}
	TAG="protein"
    elif [ ${Task} == N ] ; then
	DBTARGET=${NUCACC}
	DBDEF=${NUCDEF}
	TAG="nucleotide"
    fi
    TempDefFile=${PID}_${TAG}_TempDefDb.txt
    
    echo "Working on file "${FilePath}
    if [ -e ${FilePath}.taxo ] ; then
	rm ${FilePath}.taxo
    fi
    touch ${FilePath}.taxo
    for LINE in $(cut -f2 ${FilePath}); do
	ACC=$(echo $LINE | cut -d'|' -f4)
	echo ${ACC}
	echo ^${ACC} > ${STASKID}.${ACC}.target.txt
	while [ ! -e ${STASKID}.${ACC}.target.txt ]; do sleep 1 ; done
	
	TAXID=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${DBTARGET} | cut -f2)
	echo "^"${TAXID}"\t" > ${STASKID}.${ACC}.taxid.txt 
	while [ ! -e ${STASKID}.${ACC}.taxid.txt ]; do sleep 1 ; done
	grep -m 1 -P -f ${STASKID}.${ACC}.taxid.txt ${DBLINEAGE} > ${STASKID}.${ACC}.lineage.txt
	while [ ! -e ${STASKID}.${ACC}.lineage.txt ]; do sleep 1 ; done
	ACCorganism=$(cut -f2 ${STASKID}.${ACC}.lineage.txt)
	ACClineage=$(cut -f3 ${STASKID}.${ACC}.lineage.txt)
	ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	
	# Viruses = first element of the lineage. But for others, its "cellular organisms"
	if [ "${ACCsupKingdom}" == "cellular organisms"  ]; then
	    ACClineage=$(echo ${ACClineage} | cut -d ';' -f2-)
	    ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	fi
	
	ACCdefinition="."
	if [ "${ACCsupKingdom}" == "Viruses"  ]; then
	    #First try, on persistant file
	    ACCdefinition=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${DBDEF} | cut -f2)
	    if [ ! -z "${ACCdefinition}" ] ; then
		#Second try, on temporary file
		ACCdefinition=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${TempDefFile} | cut -f2)
		if [ ! -z "${ACCdefinition}" ] ; then        
		    echo "Unkown AccID ${ACC} in ${DBDEF}"
		    echo "Ask ebi"
		    echo "curl -s -N https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text"
		    curl -s -N "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text" > ${STASKID}.${ACC}.${TAG}.def
		    ACCdefinition=$(grep -m 1 -f DEFINITION.txt ${STASKID}.${ACC}.${TAG}.def | cut -c 13-)
		    if [ ! -z "${ACCdefinition}" ] ; then
			echo "Unable to download definition from EBI for ${ACC}"
			ACCdefinition="#N/D"
		    else
			printf "${ACC}\t${ACCdefinition}\n" >> ${TempDefFile}
		    fi
		    rm ${STASKID}.${ACC}.${TAG}.def
		fi
	    fi
	fi
	
	result="${ACC}\t${ACCorganism}\t${ACCsupKingdom}\t${ACClineage}\t${ACCdefinition}\n"
	printf "${result}" >> ${FilePath}.taxo
	
	rm ${STASKID}.${ACC}.target.txt
	rm ${STASKID}.${ACC}.taxid.txt
	rm ${STASKID}.${ACC}.lineage.txt
		
    done
done

touch Taxo_Ok/${STASKID}_Taxo.ok
