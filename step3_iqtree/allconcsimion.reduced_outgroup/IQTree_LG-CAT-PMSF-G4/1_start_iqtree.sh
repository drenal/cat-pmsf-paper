#!/bin/bash
DATASET="allconcsimion.reduced_outgroup"
DATASETSHORT="allconcsimion"
ALIGNEXT="phy"
TREEFILEEXT="fasta"

IDENTIFIER="lg_${DATASETSHORT}_icc"

iqtree2 -s ${DATASET}.${ALIGNEXT} -m LG+G4 -B 5000 --wbtl -fs ${IDENTIFIER}_chain1.sitefreq
