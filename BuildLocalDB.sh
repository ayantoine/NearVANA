#! /bin/bash

set -e #if a command crash, the script interrupt immediatly

FOLDER=$1
SELF_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
DB_SCRIPT=${SELF_SCRIPT_DIR}/Database
DIR_DATE=`date +%Y-%m-%d`

#Genbank db
if [ ! -d ${FOLDER}/${DIR_DATE} ] ; then
	mkdir ${FOLDER}/${DIR_DATE}
fi
echo "bash ${DB_SCRIPT}/DownloadData.sh ${FOLDER}/${DIR_DATE}"
bash ${DB_SCRIPT}/DownloadData.sh ${FOLDER}/${DIR_DATE}

#Viral length
echo "python ${DB_SCRIPT}/GetMeanViralLength.py -f ${FOLDER}"
python ${DB_SCRIPT}/GetMeanViralLength.py -f ${FOLDER}

#VMR file
echo "python ${DB_SCRIPT}/GetICTV-VMRfile.py -f ${FOLDER}"
python ${DB_SCRIPT}/GetICTV-VMRfile.py -f ${FOLDER}

#Check auxiliary files
if [ ! -f ${FOLDER}/NuclAccId2Def.tsv ]; then
	touch -a ${FOLDER}/NuclAccId2Def.tsv
else
	echo "${FOLDER}/NuclAccId2Def.tsv already exist. Nothing to do."
fi
if [ ! -f ${FOLDER}/ProtAccId2Def.tsv ]; then
	touch -a ${FOLDER}/ProtAccId2Def.tsv
else
	echo "${FOLDER}/ProtAccId2Def.tsv already exist. Nothing to do."
fi

