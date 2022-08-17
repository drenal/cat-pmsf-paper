#!/bin/bash
MODEL="lg"
CHAIN="bad"

for f in *.phylip; do
	BASE=`basename -s .phylip $f`

	IDENTIFIER="${MODEL}_${BASE}_${CHAIN}"

	iqtree2 -s $f -m LG+G4 -B 5000 --wbtl -fs ${IDENTIFIER}.sitefreq
done
