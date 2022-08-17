#!/bin/bash
DATASET="microsporidia"
ALIGNEXT="phy"
TREEFILEEXT="phy"

IDENTIFIER="gtr_${DATASET}_icc"

ln -s ../../../datasets/${DATASET}.${ALIGNEXT}
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.sitefreq
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.exchangeabilities

iqtree2 -s ${DATASET}.${ALIGNEXT} -m ${IDENTIFIER}_chain1.exchangeabilities+R8 -B 5000 --wbtl -fs ${IDENTIFIER}_chain1.sitefreq
