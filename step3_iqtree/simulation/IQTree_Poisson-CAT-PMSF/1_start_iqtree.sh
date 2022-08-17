#!/bin/bash
MODEL="poisson"
CHAIN="bad"
for f in *.phylip; do
	BASE=`basename -s .phylip $f`
	IDENTIFIER="${BASE}_${CHAIN}"
	iqtree2 -s $f -m Poisson+G4 -B 5000 --wbtl -fs ${IDENTIFIER}.sitefreq
done
