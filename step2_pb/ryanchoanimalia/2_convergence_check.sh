#!/bin/bash
DATASET="ryanchoanimalia"

IDENTIFIER="lg_${DATASET}_icc"
bpcomp -o ${IDENTIFIER}.bpcomp -x 3000 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2
tracecomp -x 3000 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2

IDENTIFIER="gtr_${DATASET}_icc"
bpcomp -o ${IDENTIFIER}.bpcomp -x 2700 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2
tracecomp -x 2700 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2
