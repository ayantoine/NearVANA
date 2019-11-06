#!/bin/bash

ARG=$1
source $ARG
source $CONF

TASK_ARRAY=(N X)
TargetPath=(${PID}_BlastN/${PID}_All.fa.${STASKID}.BlastN_2.tab ${PID}_BlastX/${PID}_All.fa.${STASKID}.BlastX_2.tab)

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
    touch ${FilePath}.taxo
    for ACC in $(cut -f2 ${FilePath}); do
	echo ${ACC}
	
	ACCtaxid=$(grep -m 1 -P "^[-A-Za-z0-9_.%]*\t${ACC}" ${DBTARGET} | cut -f3)
	ACCorganism=$(grep -m 1 -P "^${ACCtaxid}\t" ${DBLINEAGE} | cut -f2)
	ACClineage=$(grep -m 1 -P "^${ACCtaxid}\t" ${DBLINEAGE} | cut -f3)
	ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	# Viruses = first element of the lineage. But for others, its "cellular organisms"
	if [ "${ACCsupKingdom}" == "cellular organisms"  ]; then
	    ACClineage=$(echo ${ACClineage} | cut -d ';' -f2-)
	    ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	fi
	
	ACCdefinition="."
	if [ "${ACCsupKingdom}" == "Viruses"  ]; then
	    #First try, on persistant file
	    ACCdefinition=$(grep -m 1 -P "^${ACC}\t" ${DBDEF} | cut -f2)
	    if [ -z ${ACCdefinition} ] ; then
		#Second try, on temporary file
		ACCdefinition=$(grep -m 1 -P "^${ACC}\t" ${TempDefFile} | cut -f2)
		if [ -z ${ACCdefinition} ] ; then        
		    echo "Unkown AccID ${ACC} in ${DBLINEAGE}"
		    echo "Ask ebi"
		    echo "curl -s https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text"
		    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text" > ${Acc}.${TAG}.def
		    ACCdefinition=$(grep -m 1 "DEFINITION" ${Acc}.${TAG}.def | cut -c 13-)
		    if [ -z ${ACCdefinition} ] ; then
			echo "Unable to download definition from EBI for ${ACC}"
			ACCdefinition="#N/D"
		    else
			printf "${ACC}\t${ACCdefinition}\n" > ${TempDefFile}
		    fi
		    rm ${Acc}.${TAG}.def
		fi
	    fi
	fi
	
	result="${ACC}\t${ACCorganism}\t${ACCsupKingdom}\t${ACClineage}\t${ACCdefinition}\n"
	printf "${result}" >> ${FilePath}.taxo
		
    done
done

touch Blast${TASK}_Ok/${STASKID}_Taxo.ok
