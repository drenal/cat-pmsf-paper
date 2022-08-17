#!/bin/bash

DATASET="allconcsimion.reduced_outgroup"
DATASETSHORT="allconcsimion"
ALIGNEXT="phylip"
TREEFILEEXT="phylip"

IDENTIFIER="gtr_${DATASETSHORT}_icc"

ln -s ../../../datasets/${DATASET}.${ALIGNEXT}
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.sitefreq
ln -s ../../../step2_pb/${DATASET}/${IDENTIFIER}_chain1.exchangeabilities

