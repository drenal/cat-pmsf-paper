#!/bin/bash
DATASET="ryanchoanimalia"

IDENTIFIER="lg_${DATASET}_icc"
readpb_mpi -ss -x 3000 1 ${IDENTIFIER}_chain1
../../scripts/convert-site-dists.py ${IDENTIFIER}_chain1.siteprofiles

IDENTIFIER="gtr_${DATASET}_icc"
readpb_mpi -ss -x 2700 1 ${IDENTIFIER}_chain1
../../scripts/convert-site-dists.py ${IDENTIFIER}_chain1.siteprofiles

readpb_mpi -rr -x 2700 1 ${IDENTIFIER}_chain1
../../convert-exchangeabilities.py ${IDENTIFIER}_chain1.meanrr
