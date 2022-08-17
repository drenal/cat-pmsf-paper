#!/bin/bash
DATASET="microsporidia"
ALIGNEXT="phy"
TREEFILEEXT="phy"

IDENTIFIER="lg_${DATASET}_icc"

ln -s ../../../datasets/${DATASET}.${ALIGNEXT}
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.sitefreq

iqtree2 -s ${DATASET}.${ALIGNEXT} -m LG+R8 -B 5000 --wbtl -fs ${IDENTIFIER}_chain1.sitefreq
