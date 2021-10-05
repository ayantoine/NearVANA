#! /bin/bash

datetime1=$(date +%s)

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

ARG=$1
source $ARG
source $CONF
source $DATA
SDIR=${GITDIR}/Workflow

USE_PAIREND="$(boolean "${PAIREND}")"
USE_METADATA="$(boolean "${METADATA}")"
USE_MULTIPLEX="$(boolean "${MULTIPLEX}")"
USE_KEEPUNASSIGNED="$(boolean "${UNASSIGNED}")"

if [ ! -f ${PID}_All.Megahit_reverseAssembly.tsv ]; then
	echo -e "\t- Decompressing reverseAssembly archive"
	zcat ${PID}_All.Megahit_reverseAssembly.tsv.gz > ${PID}_All.Megahit_reverseAssembly.tsv
fi

if [ "$USE_MULTIPLEX" = true ] ; then
	if [ "$PREFILTER" = true ] ; then
		if [ "$DIAMOND" = true ] ; then
			if [ ! -f ${PID}.StatBlastD.ok ]; then
				echo "python ${SDIR}/Extract4Stat.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} > Stat_BlastD.o"
				python ${SDIR}/Extract4Stat.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} > Stat_BlastD.o
				touch ${PID}.StatBlastD.ok
			else
				echo "${PID}.StatBlastD.ok already existing, pass"
			fi
		fi
		if [ "$BLASTX" = true ] ; then
			if [ ! -f ${PID}.StatBlastX.ok ]; then
				echo "python ${SDIR}/Extract4Stat.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} > Stat_BlastX.o"
				python ${SDIR}/Extract4Stat.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} > Stat_BlastX.o
				touch ${PID}.StatBlastX.ok
			else
				echo "${PID}.StatBlastX.ok already existing, pass"
			fi
		fi
		if [ "$BLASTN" = true ] ; then
			if [ ! -f ${PID}.StatBlastN.ok ]; then
				echo "python ${SDIR}/Extract4Stat.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} > Stat_BlastN.o"
				python ${SDIR}/Extract4Stat.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} > Stat_BlastN.o
				touch ${PID}.StatBlastN.ok
			else
				echo "${PID}.StatBlastN.ok already existing, pass"
			fi
		fi
	else
		if [ "$DIAMOND" = true ] ; then
			if [ ! -f ${PID}.StatBlastD.ok ]; then
				echo "python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} -1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv > Stat_BlastD.o"
				python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv > Stat_BlastD.o
				touch ${PID}.StatBlastD.ok
			else
				echo "${PID}.StatBlastD.ok already existing, pass"
			fi
		fi
		if [ "$BLASTX" = true ] ; then
			if [ ! -f ${PID}.StatBlastX.ok ]; then
				echo "python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} -1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv > Stat_BlastX.o"
				python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv > Stat_BlastX.o
				touch ${PID}.StatBlastX.ok
			else
				echo "${PID}.StatBlastX.ok already existing, pass"
			fi
		fi
		if [ "$BLASTN" = true ] ; then
			if [ ! -f ${PID}.StatBlastN.ok ]; then
				echo "python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} -1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv > Stat_BlastN.o"
				python ${SDIR}/Extract4Stat_all.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv -d ${DATA} 1 ${LOCALDB}/All_Family_GB.list.tsv -2 ${LOCALDB}/All_Genus_GB.list.tsv -3 ${LOCALDB}/All_Species_GB.list.tsv > Stat_BlastN.o
				touch ${PID}.StatBlastN.ok
			else
				echo "${PID}.StatBlastN.ok already existing, pass"
			fi
		fi
	fi
else
	if [ "$DIAMOND" = true ] ; then
		if [ ! -f ${PID}.StatBlastD.ok ]; then
			echo "python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv > Stat_BlastD.o"
			python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastD_results.tab -o ${PID}_Stat_BlastD/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv > Stat_BlastD.o
			touch ${PID}.StatBlastD.ok
		else
			echo "${PID}.StatBlastD.ok already existing, pass"
		fi
	fi
	if [ "$BLASTX" = true ] ; then
		if [ ! -f ${PID}.StatBlastX.ok ]; then
			echo "python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv > Stat_BlastX.o"
			python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastX_results.tab -o ${PID}_Stat_BlastX/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv > Stat_BlastX.o
			touch ${PID}.StatBlastX.ok
		else
			echo "${PID}.StatBlastX.ok already existing, pass"
		fi
	fi
	if [ "$BLASTN" = true ] ; then
		if [ ! -f ${PID}.StatBlastN.ok ]; then
			echo "python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv > Stat_BlastN.o"
			python ${SDIR}/Extract4Stat_RNAseq.py -i ${PID}_BlastN_results.tab -o ${PID}_Stat_BlastN/ -v ${VMR} -r ${PID}_All.Megahit_reverseAssembly.tsv > Stat_BlastN.o
			touch ${PID}.StatBlastN.ok
		else
			echo "${PID}.StatBlastN.ok already existing, pass"
		fi
	fi
fi

if [ -f ${PID}_All.Megahit_reverseAssembly.tsv ]; then
	echo -e "\t- remove reverseAssembly file"
	rm ${PID}_All.Megahit_reverseAssembly.tsv
fi


touch ${PID}.BasicStat.ok

datetime2=$(date +%s)
delta=$((datetime2 - datetime1))
echo "Time BasicStat: "$delta > Time04.txt
