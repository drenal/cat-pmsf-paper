#!/bin/bash
../scripts/filter.py -i allconcsimion.fasta -b blacklist.txt -o allconcsimion.reduced_outgroup.fasta
../scripts/split_msa.py -i allconcsimion.reduced_outgroup.fasta -t fasta -s phylip -o allconcsimion.reduced_outgroup -p 1
mv allconcsimion.reduced_outgroup_0.phylip allconcsimion.reduced_outgroup.phylip
