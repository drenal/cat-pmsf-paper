#!/bin/bash
DATASET="ryanchoanimalia"
ALIGNEXT="phy"
TREEFILEEXT="fasta"

IDENTIFIER="gtr_${DATASET}_icc"

iqtree2 -s ${DATASET}.${ALIGNEXT} -m ${IDENTIFIER}_chain1.exchangeabilities+G4 -B 5000 --wbtl -fs ${IDENTIFIER}_chain1.sitefreq
