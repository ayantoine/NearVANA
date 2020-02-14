#!/bin/bash

ARG=$1
source $ARG
source $CONF



FilePath=${PID}_Diamond/${PID}_All.fa.${STASKID}.Diamond_2.tab

TempDefFile=${PID}_protein_TempDefDb.txt

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
    
    TAXID=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${PROACC} | cut -f2)
    echo "^"${TAXID}"\t" > ${STASKID}.${ACC}.taxid.txt 
    while [ ! -e ${STASKID}.${ACC}.taxid.txt ]; do sleep 1 ; done
    grep -m 1 -P -f ${STASKID}.${ACC}.taxid.txt ${DBLINEAGE} > ${STASKID}.${ACC}.lineage.txt
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
	ACCdefinition=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${PRODEF} | cut -f2)
	ACCsize=${#ACCdefinition}
	if [ "$ACCsize" -le 1 ] ; then #less than or equal to
	    #Second try, on temporary file
	    ACCdefinition=$(grep -m 1 -f ${STASKID}.${ACC}.target.txt ${TempDefFile} | cut -f2)
	    ACCsize=${#ACCdefinition}
	    if [ "$ACCsize" -le 1 ] ; then #less than or equal to
		echo "Unkown AccID ${ACC} in ${PRODEF}"
		echo "Ask ebi"
		echo "curl -s -N https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ACC}&rettype=gb&retmode=text"
		curl -s -N "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ACC}&rettype=gb&retmode=text" > ${STASKID}.${ACC}.protein.def
		ACCdefinition=$(grep -A 1 -m 1 -f DEFINITION.txt ${STASKID}.${ACC}.protein.def | grep -v -f ACCESSION.txt | cut -c 13- | tr '\n' ' ')
		ACCsize=${#ACCdefinition}
		if [ "$ACCsize" -le 1 ] ; then #less than or equal to
		    echo "Unable to download definition from EBI for ${ACC}"
		    ACCdefinition="#N/D"
		else
		    printf "${ACC}\t${ACCdefinition}\n" >> ${TempDefFile}
		fi
		rm ${STASKID}.${ACC}.protein.def
	    fi
	fi
    fi
    
    result="${ACC}\t${ACCorganism}\t${ACCsupKingdom}\t${ACClineage}\t${ACCdefinition}\n"
    printf "${result}" >> ${FilePath}.taxo
    
    rm ${STASKID}.${ACC}.target.txt
    rm ${STASKID}.${ACC}.taxid.txt
    rm ${STASKID}.${ACC}.lineage.txt
	    
done


touch Taxo_Ok/${STASKID}_Taxo.ok
