#!/bin/bash
DATASET="allconcsimion.reduced_outgroup"

ln -s ../../datasets/${DATASET}.phylip
ln -s ../../step1_iqtree_lg/${DATASET}/${DATASET}.*.treefile
