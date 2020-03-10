#!/bin/bash

ARG=$1
source $ARG
source $CONF

nb_jobs=$(ls ${PID}_ToBlast | wc -l)

if [ ! -f Diamond_Ok/${STASKID}_Diamond.ok ]; then
    echo "------ Viral database blast ------"
    echo "diamond blastx -d ${VIRPTDB_DIAMOND} -q ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -o ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.Diamond_1.tab -p ${MULTICPU} --max-target-seqs 1 --max-hsps 1 --evalue 0.001 --block-size $((${MULTIMEMORY}/6))"
    diamond blastx -d ${VIRPTDB_DIAMOND} -q ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -o ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.Diamond_1.tab -p ${MULTICPU} --max-target-seqs 1 --max-hsps 1 --evalue 0.001 --block-size $((${MULTIMEMORY}/6))
    python ${SDIR}/SplitFastaOnBlastResults.py -f ${PID}"_ToBlast"/${PID}_All.fa.${STASKID} -b ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.Diamond_1.tab -o ${PID}"_BlastD"/${PID}_All.fa.${STASKID}
    echo "------ /Viral database blast ------"

    echo "------ Complete database blast ------"
    echo "diamond blastx -d ${ALLPTDB_DIAMOND} -q ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.keeped -o ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.Diamond_2.tab -p ${MULTICPU} --max-target-seqs 5 --max-hsps 1 --evalue 0.001 --block-size $((${MULTIMEMORY}/6))"
    diamond blastx -d ${ALLPTDB_DIAMOND} -q ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.keeped -o ${PID}"_BlastD"/${PID}_All.fa.${STASKID}.Diamond_2.tab -p ${MULTICPU} --max-target-seqs 5 --max-hsps 1 --evalue 0.001 --block-size $((${MULTIMEMORY}/6))
    echo "------ /Complete database blast ------"

    touch Diamond_Ok/${STASKID}_Diamond.ok
else
    echo "File Diamond_Ok/Diamond${TASK}.ok already present. Assume there is no need to launch Diamond on this index-array."
fi
