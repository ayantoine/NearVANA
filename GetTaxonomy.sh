#!/bin/bash

datetime1=$(date +%s)

ARG=$1
source $ARG
source $CONF

Task=$2

if [ ! -f Taxo${Task}_Ok/${STASKID}_Taxo.ok ] ; then
    if [ ${Task} == X ] ; then
	DBTARGET=${PROACC}
	DBDEF=${PRODEF}
	TAG="protein"
	FilePath=${PID}_Blast${Task}/${PID}_All.fa.${STASKID}.Blast${Task}_2.tab
    elif [ ${Task} == N ] ; then
	DBTARGET=${NUCACC}
	DBDEF=${NUCDEF}
	TAG="nucleotide"
	FilePath=${PID}_Blast${Task}/${PID}_All.fa.${STASKID}.Blast${Task}_2.tab
    elif [ ${Task} == D ] ; then
	DBTARGET=${PROACC}
	DBDEF=${PRODEF}
	TAG="protein"
	FilePath=${PID}_Blast${Task}/${PID}_All.fa.${STASKID}.Diamond_2.tab
    fi
    TempDefFile=${PID}_${Task}_TempDefDb.txt
    
    echo "Working on file "${FilePath}
    if [ -e ${FilePath}.taxo ] ; then
	rm ${FilePath}.taxo
    fi
    touch ${FilePath}.taxo
    for LINE in $(cut -f2 ${FilePath} | sort -u); do
	ACC=$(echo $LINE | cut -d'|' -f4)
	echo ${ACC}
	echo ^${ACC} > ${STASKID}.${ACC}.target.txt
	while [ ! -e ${STASKID}.${ACC}.target.txt ]; do sleep 1 ; done
	
	TAXID=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${DBTARGET} | cut -f2)
	echo ${ACC}"\t"${TAXID}
	echo "^"${TAXID}"\t" > ${STASKID}.${ACC}.taxid.txt 
	while [ ! -e ${STASKID}.${ACC}.taxid.txt ]; do sleep 1 ; done
	LINEAGE=$(grep -m 1 -P -f ${STASKID}.${ACC}.taxid.txt ${DBLINEAGE})
	echo ${ACC}"\t"${LINEAGE}
	echo ${LINEAGE} > ${STASKID}.${ACC}.lineage.txt
	while [ ! -e ${STASKID}.${ACC}.lineage.txt ]; do sleep 1 ; done
	ACCorganism=$(cut -f2 ${STASKID}.${ACC}.lineage.txt)
	ACClineage=$(cut -f3 ${STASKID}.${ACC}.lineage.txt)
	ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	
	# Viruses = first element of the lineage. But for others, its "cellular organisms"
	if [ "${ACCsupKingdom}" == "cellular organisms"  ]; then
	    ACClineage=$(echo ${ACClineage} | cut -d ';' -f2- | sed -e 's/^ //')
	    ACCsupKingdom=$(echo ${ACClineage} | cut -d ';' -f1)
	fi
	
	ACCdefinition="."
	if [ "${ACCsupKingdom}" == "Viruses"  ]; then
	    #First try, on persistant file
	    ACCdefinition=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${DBDEF} | cut -f2)
	    ACCsize=${#ACCdefinition}
	    if [ "$ACCsize" -le 1 ] ; then #less than or equal to
		#Second try, on temporary file
		ACCdefinition=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${TempDefFile} | cut -f2)
		ACCsize=${#ACCdefinition}
		if [ "$ACCsize" -le 1 ] ; then #less than or equal to
		    echo "Unkown AccID ${ACC} in ${DBDEF}"
		    echo "Ask ebi"
		    echo "curl -s -N https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text"
		    curl -s -N "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${TAG}&id=${ACC}&rettype=gb&retmode=text" > ${STASKID}.${ACC}.${TAG}.def
		    ACCdefinition=$(grep -A 1 -m 1 -f DEFINITION.txt ${STASKID}.${ACC}.${TAG}.def | grep -v -f ACCESSION.txt | cut -c 13- | tr '\n' ' ')
		    ACCsize=${#ACCdefinition}
		    if [ "$ACCsize" -le 1 ] ; then #less than or equal to
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
    
    touch Taxo${Task}_Ok/${STASKID}_Taxo.ok
else
    echo "Taxo${Task}_Ok/${STASKID}_Taxo.ok already existing, do nothing..."
fi

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time Taxo${Task}: "$delta #> Time11-${Task}.txt
