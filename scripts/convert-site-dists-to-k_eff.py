#!/usr/bin/env python3
"""
Convert IQTree's PMSF site distributions 
to effective number of amino acids.

Authors:
- Dominik Schrempf
- Lenard Szantho
"""

import argparse
import os
import numpy as np
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('infile', help="IQ-TREE site distribution file.")
    
    return parser.parse_args()


def np_array_to_string(a, d=' '):
    """Convert a numpy array to a string."""
    x_arrstr = np.char.mod('%.8f', a)
    return d.join(x_arrstr)


def p_homo(a):
    """Homoplasy."""
    # return sum([x**2 for x in a])
    return np.sum(np.square(a))

def p_entropy(a):
    """ Entropy. """
    return (-1) * np.sum(a * np.log(a))

def keff_entropy(a):
    """ Effective number of amino acids (entropy). """
    return np.exp(p_entropy(a))

def keff_homoplasy(a):
    """Effective number of amino acids (homoplasy)."""
    return 1.0/np.sum(np.square((a)))


a = parse_arguments()

fn = a.infile

of = os.path.splitext(fn)[0] + ".keff"

df = pd.read_table(fn, sep=r'\s+', header=None, index_col=0)

df.index.names = ['Site']

df_entropy = df.apply(func=keff_entropy, axis=1)
df_entropy.to_csv(of+"Entropy")

df_homoplasy = df.apply(func=keff_homoplasy, axis=1)
df_homoplasy.to_csv(of)

