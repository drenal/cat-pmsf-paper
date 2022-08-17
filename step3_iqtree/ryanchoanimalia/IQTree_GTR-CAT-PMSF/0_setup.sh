#!/bin/bash

DATASET="ryanchoanimalia"
ALIGNEXT="phy"
TREEFILEEXT="phy"

IDENTIFIER="gtr_${DATASET}_icc"

ln -s ../../../datasets/${DATASET}.${ALIGNEXT}
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.sitefreq
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.exchangeabilities

