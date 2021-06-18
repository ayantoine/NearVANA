#! /bin/bash

set -e #if a command crash, the script interrupt immediatly

FOLDER=$1
SELF_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
DB_SCRIPT=${SELF_SCRIPT_DIR}/Database
DIR_DATE=`date +%Y-%m-%d`

#Genbank db
mkdir ${FOLDER}/${DIR_DATE}
echo "bash ${DB_SCRIPT}/DownloadData.sh ${FOLDER}/${DIR_DATE}"
bash ${DB_SCRIPT}/DownloadData.sh ${FOLDER}/${DIR_DATE}

#Viral length
echo "python ${DB_SCRIPT}/GetMeanViralLength.py -f ${FOLDER}"
python ${DB_SCRIPT}/GetMeanViralLength.py -f ${FOLDER}

#Check auxiliary files
touch -a ${FOLDER}/NuclAccId2Def.tsv
touch -a ${FOLDER}/ProtAccId2Def.tsv
