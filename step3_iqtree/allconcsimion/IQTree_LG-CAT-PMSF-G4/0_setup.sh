#!/bin/bash

DATASET="allconcsimion"
ALIGNEXT="phy"
TREEFILEEXT="fasta"

IDENTIFIER="lg_${DATASET}_icc"

ln -s ../../../datasets/${DATASET}.${ALIGNEXT}
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.sitefreq

