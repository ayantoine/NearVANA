#! /bin/bash

set -e #if a command crash, the script interrupt immediatly

FOLDER=$1

FTPlineage="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/"
TarLineage="new_taxdump"
SubFileLineage="fullnamelineage.dmp"
SubFileNodes="nodes.dmp"
SubFileNames="names.dmp"

FTPaccId="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid"
NTarget="nucl_gb.accession2taxid"
PTarget="prot.accession2taxid"
TargetArray=(${NTarget} ${PTarget})

SELF_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

##Lineage part
#Download
echo ${FTPlineage}/${TarLineage}.tar.gz.md5 --output ${FOLDER}/${TarLineage}.md5
curl ${FTPlineage}/${TarLineage}.tar.gz.md5 --output ${FOLDER}/${TarLineage}.md5
echo ${FTPlineage}/${TarLineage}.tar.gz --output ${FOLDER}/${TarLineage}.tar.gz
curl ${FTPlineage}/${TarLineage}.tar.gz --output ${FOLDER}/${TarLineage}.tar.gz

#Checksum verification
echo "------CheckSum verification------"
CheckRef=`cut -f1 -d" " ${FOLDER}/${TarLineage}.md5`
echo "CheckRef ${FOLDER}/${TarLineage}.md5 $CheckRef"
md5sum ${FOLDER}/${TarLineage}.tar.gz > ${FOLDER}/${TarLineage}.current.md5
CheckSum=`cut -f1 -d" " ${FOLDER}/${TarLineage}.current.md5`
echo "CheckSum ${FOLDER}/${TarLineage}.current.md5 $CheckSum"
echo "------/CheckSum verification------"

if test "$CheckRef" != "$CheckSum" ; then
	#rm ${FOLDER}/${TarLineage}.tar.gz
	echo "Unable to dowload accurate file for "${FOLDER}/${TarLineage}
	exit 1
fi

#extract specific file
tar -xf ${FOLDER}/${TarLineage}.tar.gz ${SubFileLineage}
mv ${SubFileLineage} ${FOLDER}/${SubFileLineage}
tar -xf ${FOLDER}/${TarLineage}.tar.gz ${SubFileNodes}
mv ${SubFileNodes} ${FOLDER}/${SubFileNodes}
tar -xf ${FOLDER}/${TarLineage}.tar.gz ${SubFileNames}
mv ${SubFileNames} ${FOLDER}/${SubFileNames}

#Format lineage file
scp ${FOLDER}/${SubFileLineage} ${FOLDER}/TempFile.txt
cut -f 1,3,5 ${FOLDER}/TempFile.txt > ${FOLDER}/${SubFileLineage}
rm ${FOLDER}/TempFile.txt

#Extract data
python ${SELF_SCRIPT_DIR}/ExtractTaxonomyFromGenbankData.py -1 ${FOLDER}/${SubFileNodes} -2 ${FOLDER}/${SubFileNames} -t family -o ${FOLDER}/All_Family_GB.list.tsv
python ${SELF_SCRIPT_DIR}/ExtractTaxonomyFromGenbankData.py -1 ${FOLDER}/${SubFileNodes} -2 ${FOLDER}/${SubFileNames} -t genus -o ${FOLDER}/All_Genus_GB.list.tsv
python ${SELF_SCRIPT_DIR}/ExtractTaxonomyFromGenbankData.py -1 ${FOLDER}/${SubFileNodes} -2 ${FOLDER}/${SubFileNames} -t species -o ${FOLDER}/All_Species_GB.list.tsv

#Store new checksum, then remove file
rm ${FOLDER}/${SubFileNodes}
rm ${FOLDER}/${SubFileNames}
rm ${FOLDER}/${TarLineage}.tar.gz
rm ${FOLDER}/${TarLineage}.md5
rm ${FOLDER}/${TarLineage}.current.md5


for Target in "${TargetArray[@]}"; do
	#Download checksum file
	echo ${FTPaccId}/${Target}.gz.md5 --output ${FOLDER}/${Target}.md5
    curl ${FTPaccId}/${Target}.gz.md5 --output ${FOLDER}/${Target}.md5
    echo ${FTPaccId}/${Target}.gz --output ${FOLDER}/${Target}.gz
    curl ${FTPaccId}/${Target}.gz --output ${FOLDER}/${Target}.gz
	
	#Checksum verification
	echo "------CheckSum verification------"
	CheckRef=`cut -f1 -d" " ${FOLDER}/${Target}.md5`
	echo "CheckRef ${FOLDER}/${Target}.md5 $CheckRef"
	md5sum ${FOLDER}/${Target}.gz > ${FOLDER}/${Target}.current.md5
	CheckSum=`cut -f1 -d" " ${FOLDER}/${Target}.current.md5`
	echo "CheckSum ${FOLDER}/${Target}.current.md5 $CheckSum"
	echo "------/CheckSum verification------"
	
	if test "$CheckRef" != "$CheckSum" ; then
		#rm ${FOLDER}/${Target}.gz
		echo "Unable to dowload accurate file for "${FOLDER}/${Target}
		exit 1
	fi
	
	#extract file
	gunzip ${FOLDER}/${Target}.gz
	
	#Remove checksum file
	rm ${FOLDER}/${Target}.current.md5
	rm ${FOLDER}/${Target}.md5
	
	#Create short version
	cut -f2,3 ${FOLDER}/${Target} > ${FOLDER}/Viral_${Target}.short.tsv
	rm ${FOLDER}/${Target}

done

##Request by: grep -m 1 "NZ_CP022958" nucl_gb.accession2taxid | cut -f 3 | grep -w -f /dev/stdin fullnamelineage.dmp 
#if [ $2 == "N" ]
	#then
	#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=gb&retmode=txt" > Temp4_$1.tab
#fi
#if [ $2 == "P" ]
	#then
	#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ACC}&rettype=gb&retmode=txt" > Temp4_$1.tab
#fi
# grep -m 1 "DEFINITION" Test.txt | cut -c 13-
# echo $VARTAXO | cut -d ';' -f1

