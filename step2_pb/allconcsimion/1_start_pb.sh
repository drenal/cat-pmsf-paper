#!/bin/bash

DATASET="allconcsimion"
ALIGNEXT="phy"
TREEFILEEXT="fasta"

for chain in {1,2}; do
    for model in {lg,gtr}; do
        mpirun -np 48 pb_mpi -cat -${model} -d ${DATASET}.${ALIGNEXT} -T ${DATASET}.${TREEFILEEXT}.treefile ${model}_${DATASET}_icc_chain${chain}
    done
done
