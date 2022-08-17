#!/bin/bash

DATASET="simulation/felsenstein*"
ALIGNEXT="phylip"
TREEFILEEXT="phylip"

model="gtr"

for f in ../../../datasets/simulation/*.${ALIGNEXT}; do
	BASE=`basename -s .phylip $f`
	ln -s $f

	chain="bad"
	IDENTIFIER="${model}_${BASE}_${chain}"
	ln -s ../../../step2_pb/simulation/${IDENTIFIER}.sitefreq
	if [[ "$model" == "gtr" ]]; then
		ln -s ../../../step2_pb/simulation/${IDENTIFIER}.exchangeabilities
	fi
	if [[ "$model" == "" ]]; then
		chain="good"
		IDENTIFIER="${model}_${BASE}_${chain}"
		ln -s ../../../step2_pb/simulation/${IDENTIFIER}.sitefreq
	fi
done

