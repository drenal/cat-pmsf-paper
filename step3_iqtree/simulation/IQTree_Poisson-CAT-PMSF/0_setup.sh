#!/bin/bash

model=""

for f in ../../../datasets/simulation/felsenstein*.phylip; do
	BASE=`basename -s .phylip $f`
	ln -s $f

	chain="bad"
	IDENTIFIER="${model}${BASE}_${chain}"
	ln -s ../../../step2_pb/simulation/${IDENTIFIER}.sitefreq
	if [[ "$model" == "gtr" ]]; then
		ln -s ../../../step2_pb/simulation/${IDENTIFIER}.exchangeabilities
	fi
done

