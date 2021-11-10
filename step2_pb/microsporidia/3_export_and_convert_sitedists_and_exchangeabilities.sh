#!/bin/bash
DATASET="microsporidia"

IDENTIFIER="lg_${DATASET}_icc"
readpb_mpi -ss -x 4000 1 ${IDENTIFIER}_chain1
../../scripts/convert-site-dists.py ${IDENTIFIER}_chain1.siteprofiles

IDENTIFIER="gtr_${DATASET}_icc"
readpb_mpi -ss -x 1000 1 ${IDENTIFIER}_chain1
../../scripts/convert-site-dists.py ${IDENTIFIER}_chain1.siteprofiles
    
readpb_mpi -rr -x 1000 1 ${IDENTIFIER}_chain1
../../convert-exchangeabilities.py ${IDENTIFIER}_chain1.meanrr
