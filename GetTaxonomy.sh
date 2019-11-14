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
	echo ^${ACC} > ${ACC}.target.txt
	while [ ! -e ${ACC}.target.txt ]; do sleep 1 ; done
	
	TAXID=$(grep -m 1 -f ${ACC}.target.txt ${DBTARGET} | cut -f2)
	echo "^"${TAXID}"\t" > ${ACC}.taxid.txt 
	grep -m 1 -P -f ${ACC}.taxid.txt ${DBLINEAGE} > ${ACC}.lineage.txt
	while [ ! -e ${ACC}.lineage.txt ]; do sleep 1 ; done
	ACCorganism=$(cut -f2 ${ACC}.lineage.txt)
	ACClineage=$(cut -f3 ${ACC}.lineage.txt)
	ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	
	# Viruses = first element of the lineage. But for others, its "cellular organisms"
	if [ "${ACCsupKingdom}" == "cellular organisms"  ]; then
	    ACClineage=$(echo ${ACClineage} | cut -d ';' -f2-)
	    ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	fi
	
	ACCdefinition="."
	if [ "${ACCsupKingdom}" == "Viruses"  ]; then
	    #First try, on persistant file
	    ACCdefinition=$(grep -m 1 -f ${ACC}.target.txt ${DBDEF} | cut -f2)
	    if [ ! -n ${ACCdefinition} ] ; then
		#Second try, on temporary file
		ACCdefinition=$(grep -m 1 -f ${ACC}.target.txt ${TempDefFile} | cut -f2)
		if [ ! -n ${ACCdefinition} ] ; then        
		    echo "Unkown AccID ${ACC} in ${DBDEF}"
		    echo "Ask ebi"
		    echo "curl -s -N https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text"
		    curl -s -N "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text" > ${ACC}.${TAG}.def
		    ACCdefinition=$(grep -m 1 -f DEFINITION.txt ${ACC}.${TAG}.def | cut -c 13-)
		    if [ ! -n ${ACCdefinition} ] ; then
			echo "Unable to download definition from EBI for ${ACC}"
			ACCdefinition="#N/D"
		    else
			printf "${ACC}\t${ACCdefinition}\n" >> ${TempDefFile}
		    fi
		    rm ${ACC}.${TAG}.def
		fi
	    fi
	fi
	
	result="${ACC}\t${ACCorganism}\t${ACCsupKingdom}\t${ACClineage}\t${ACCdefinition}\n"
	printf "${result}" >> ${FilePath}.taxo
	
	rm ${ACC}.target.txt
	rm ${ACC}.taxid.txt
	rm ${ACC}.lineage.txt
		
    done
done

touch Taxo_Ok/${STASKID}_Taxo.ok
