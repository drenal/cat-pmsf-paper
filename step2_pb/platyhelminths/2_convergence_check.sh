#!/bin/bash
DATASET="platyhelminths"

IDENTIFIER="lg_${DATASET}_icc"
bpcomp -o ${IDENTIFIER}.bpcomp -x 1000 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2
tracecomp -x 1000 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2

IDENTIFIER="gtr_${DATASET}_icc"
bpcomp -o ${IDENTIFIER}.bpcomp -x 1000 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2
tracecomp -x 1000 ${IDENTIFIER}_chain1 ${IDENTIFIER}_chain2
