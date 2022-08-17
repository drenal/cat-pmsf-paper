#!/bin/bash
DATASET="ryanchoanimalia"
ALIGNEXT="phy"
TREEFILEEXT="fasta"

IDENTIFIER="lg_${DATASET}_icc"

iqtree2 -s ${DATASET}.${ALIGNEXT} -m LG+G4 -B 5000 --wbtl -fs ${IDENTIFIER}_chain1.sitefreq
