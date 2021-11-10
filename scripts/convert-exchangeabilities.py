#!/usr/bin/env python3
"""
This script transforms the amino acid exchangeabilities 
infered by Phylobayes and outputted by readpb_mpi utility
with the -rr switch into PAML format exchangeability 
matrix which is expected by IQTree if one wants to specify
a custom exchangeability model.
"""

import argparse

import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('exchangeabilities', help='Path to the exchangeabilities file created by readpb_mpi (file extensions must end with ".meanrr")')
    args = parser.parse_args()

    output_file_name = args.exchangeabilities.rstrip('.meanrr') + ".exchangeabilities"

    exchangeabilities = read_phylobayes_exchangeabilities(args.exchangeabilities)

    write_PAML_exchangeabilities(output_file_name, exchangeabilities)


def read_phylobayes_exchangeabilities(filename):
    """
Parse a readpb_mpi created exchangeability file.

Phylobayes outputs a file with the first line specifying the
order of the amino acids, e.g.:

A C D E F G H I K L M N P Q R S T V W Y 

then an empty line and then 3 columns lines start
listing the exchangeabilities:

A	C	1.37435
A	D	0.486164
A	E	1.23747
A	F	0.20664
...

Parameters
----------
filename : string
    Path to .meanrr file

Returns
-------
exchangeabilities : dict
    Dictionary of dictionary contianing the exchangeabilities
    between pairs of amino acids.
    {"A": {"C": 0.5, "D": 0.3, ...}, "C": {"D": 0.3, ...}}
oder_of_AA : array
    Order of amino acids in which they were listed in the file
    
    """
    
    with open(filename, mode='r') as fh:
        order_of_AA = fh.readline().split()
        exchangeabilities = pd.DataFrame(
            np.zeros((len(order_of_AA), len(order_of_AA)), dtype=np.float64),
            index=order_of_AA,
            columns=order_of_AA,
            dtype=np.float64)
        
        for line in fh:
            if len(line.strip()) > 0:
                source, target, exch = line.split()
                exchangeabilities[source][target] = float(exch)
                exchangeabilities[target][source] = float(exch)

    return exchangeabilities

        
def write_PAML_exchangeabilities(filename, exchangeabilities):
    """
Write exchangeabilities to filename in PAML (Phylogenetic Analysis by 
Maximum Likelihood software package written by Ziheng Yang) format.

PAML format file contains the exchangeability matrix in a strict way:
- the order of elements is: 
A    R    N    D    C    Q    E   G   H    I    L    K    M    F    P    S    T    W    Y    V
Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val
- only the lower diagonal triangle of the matrix is listed

Parameters
----------
filename: string
    Path to output file to write to
exchangeabilities : dict
    Dictionary of dictionary contianing the exchangeabilities
    between pairs of amino acids.
    {"A": {"C": 0.5, "D": 0.3, ...}, "C": {"D": 0.3, ...}}
    """
    
    PAML_order_of_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S' , 'T', 'W', 'Y', 'V']

    # calculate the diagonal elements
    # according to Lartillot 2006
    # Q_{aa} = sum_{b \neq a} Q_{ab}
    
    # naive solution
    for col in exchangeabilities.columns:
        exchangeabilities[col][col] = -1.0 * exchangeabilities.loc[:,col].sum()

    # TODO:
    # one can set value for an entire column:
    # df.loc[:, 'max_speed'] = 30
    # so if we can provide a dataframe in the right order on the right side,
    # it could be made nicer

    # TODO:
    # make this write nicer...
    with open(filename, mode='w') as fh:
        for source_index, source_aa in enumerate(PAML_order_of_AA):
            for target_aa in PAML_order_of_AA[:source_index]:
                fh.write("{:08.6f} ".format(exchangeabilities[source_aa][target_aa]))
            fh.write('\n')
        fh.write('\n' + '0.050000 '*20 + '\n')


if __name__ == '__main__':
    main()


