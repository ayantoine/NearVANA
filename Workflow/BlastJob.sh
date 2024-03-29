#!/bin/bash

ARG=$1
source $ARG
source $CONF
SDIR=${GITDIR}/Workflow

function boolean() {
  case $1 in
    TRUE) echo true ;;
    FALSE) echo false ;;
    *) echo "Err: Unknown boolean value \"$1\"" 1>&2; exit 1 ;;
   esac
}

PREFILTER="$(boolean "${PREFILTER}")"

#N or X
TASK=$2

if [ ${TASK} == N ]; then
    BLAST="blastn"
    BLAST_OPT="-dust yes -task blastn"
    VIRDB=${VIRNTDB}
    ALLDB=${ALLNTDB}
    TAG="nucleotide"
elif [ ${TASK} == X ]; then
    BLAST="blastx"
    BLAST_OPT="-seg yes -matrix BLOSUM62"
    VIRDB=${VIRPTDB}
    ALLDB=${ALLPTDB}
    TAG="protein"
fi

nb_jobs=$(ls ${PID}_ToBlast | wc -l)

if [ ! -f Blast${TASK}_Ok/${STASKID}_Blast${TASK}.ok ]; then
    if [ "$PREFILTER" = true ] ; then
        echo "------ Viral database blast ------"
        echo "${BLAST} ${BLAST_OPT} -strand both -query ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -db ${VIRDB} -evalue 0.001 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_1.tab"
        ${BLAST} ${BLAST_OPT} -strand both -query ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -db ${VIRDB} -evalue 0.001 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_1.tab
        python ${SDIR}/SplitFastaOnBlastResults.py -f ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -b ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_1.tab -o ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}
        echo "------ /Viral database blast ------"

        echo "------ Complete database blast ------"
        echo "${BLAST} ${BLAST_OPT} -strand both -query ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.keeped -db ${ALLDB} -evalue 0.001 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_2.tab"
        ${BLAST} ${BLAST_OPT} -strand both -query ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.keeped -db ${ALLDB} -evalue 0.001 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_2.tab
        echo "------ /Complete database blast ------"
    else
        echo "------ Complete database blast ------"
        touch ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_1.tab
        scp ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.keeped
        echo "${BLAST} ${BLAST_OPT} -strand both -query ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -db ${ALLDB} -evalue 0.001 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_2.tab"
        ${BLAST} ${BLAST_OPT} -strand both -query ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -db ${ALLDB} -evalue 0.001 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out ${PID}"_Blast${TASK}"/${PID}_All.fa.${STASKID}.Blast${TASK}_2.tab
        echo "------ /Complete database blast ------"
    fi

    touch Blast${TASK}_Ok/${STASKID}_Blast${TASK}.ok
else
    echo "File Blast${TASK}_Ok/${STASKID}_Blast${TASK}.ok already present. Assume there is no need to launch Blast on this index-array."
fi
