#!/usr/bin/env python3

import argparse

import os

import numpy as np

d = """Convert PhyloBayes site distributions to IQ-TREE site frequencies.

1. Obtain PhyloBayes site distributions from a CAT run with chain name

`CHAIN_NAME`, and a burn in of 100 MCMC jumps:

    readpb_mpi -ss -x 100 1 CHAIN_NAME

2. Convert the site distributions `CHAIN_NAME.siteprofiles` to IQ-TREE format

with this script:

    convert-site-dists CHAIN_NAME.siteprofiles

3. Use the site frequency file "CHAIN_NAME.sitefreq" with IQ-TREE and PMSF:

    iqtree -s ALIGNMENT -m LG+R4 -fs CHAIN_NAME.sitefreq     (? DOES NOT WORK)

    iqtree -s ALIGNMENT -m LG+C10+R4 -fs CHAIN_NAME.sitefreq (THIS SEEMS TO WORK, but C10?)

"""

def parse_arguments():

    parser = argparse.ArgumentParser(description=d)

    parser.add_argument('infile', help="PhyloBayes site distribution file.")

    return parser.parse_args()

a = parse_arguments()

fn = a.infile

of = os.path.splitext(fn)[0] + ".sitefreq"

def pb_to_iqtree(pb):

    amino_acids_pb = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",

                      "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    amino_acids_iqtree = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',

                          'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    amino_acids_pb_dict = dict((a, i) for i, a in enumerate(amino_acids_pb))

    amino_acids_iqtree_indices = [amino_acids_pb_dict[a] for a in

                                  amino_acids_iqtree]

    iqtree = np.array([pb[i] for i in amino_acids_iqtree_indices])

    return iqtree

def np_array_to_string(a, d=' '):

    """Convert a numpy array to a string."""

    x_arrstr = np.char.mod('%.8f', a)

    return d.join(x_arrstr)

def repair(a):

    """Remove zeroes."""

    eps = 1e-8

    return np.array([x if (x >= eps) else eps for x in a])

def normalize(a):

    """Normalize stationary distributions."""

    return a / a.sum()

ofh = open(of, mode='x')

with open(fn, mode='r') as fo:

    # Skip two header lines.

    fo.readline()

    fo.readline()

    for ln in fo:

        fields = ln.strip().split()

        pos = int(fields[0])

        dist = np.array([float(e) for e in fields[1:]])

        iqtree = pb_to_iqtree(dist)

        iqtree_norm = normalize(repair(iqtree))

        print(pos, np_array_to_string(iqtree_norm), file=ofh)

ofh.close()
