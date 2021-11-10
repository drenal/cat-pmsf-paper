#!/usr/bin/env python3
"""
This script removes the sequences specified in a blocklist file by
sequence IDs from a fasta format file. Now it also The blocklist file can contain the same ID formats as the fasta file or only the gisaid ID number. Use the appropiate positional arguments to choose between this two different cases.

Authors:
- Lenard Szantho <lenard@drenal.eu>
- Mario Perez Jimenez <marioperez@ciencias.unam.mx>

Version: 0.2 (2020-10-27)
"""

import argparse
import sys
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        action="store",
        help="Path to fasta file",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--blocklist",
        dest="blocklist",
        action="store",
        help="Path to blocklist file containing sequence IDs to be removed from the specified fasta file",
        required=False,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        action="store",
        help="Path to the output file",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--allowlist",
        dest="allowlist",
        action="store_true",
        help="Make allowlist from blocklist. Instead of filtering out blocklist ones, only print out their sequences",
        default=False,
    )
    parser.add_argument(
        "-d",
        "--delete",
        dest="delete",
        action="store",
        help="Replace these chars by gap",
        default="",
        required=False,
    )

    parser.add_argument(
        "-sep",
        "--sep",
        dest="separator",
        action="store",
        help="Separator for split the ID's in the fasta file",
        required=False,
    )

    parser.add_argument(
        "-g",
        "--gisaid_id",
        dest="gis_id",
        action="store_true",
        help="Use it when the blocklist file contains only gisaid ids",
        required=False,
    )

    args = parser.parse_args()

    sequences = pd.read_csv(args.input, lineterminator=">", names=["Seq"])

    sequences = pd.DataFrame(
        sequences["Seq"].str.split(r"\n|\r", 1).tolist(), columns=["ID", "Sequence"]
    )
    if args.gis_id:
        if args.separator == None:
            print("This approach needs a separator to split the fasta file ID's")
            sys.exit(0)
        else:
            metadata = pd.DataFrame(
                sequences["ID"].str.split(args.separator, 3).tolist(), columns = ["gisaid_id","date","location"]
            )

            sequences = sequences.merge(metadata, left_index=True, right_index=True)

    def filter_function(args_blocklist,dataf,args_gis_id,args_allow):
        if args_gis_id:
            col_id = 'gisaid_id'
        else:
            col_id = 'ID'
        if args_allow:
            blocklist_num = len(dataf[dataf[col_id].isin(blocklist["ID"])])
            if blocklist_num == 0:
                print('There is no mathcing sequences to allow, check that the format of the blocklist and fasta file have the same format')
            dataf = dataf[dataf[col_id].isin(blocklist["ID"])]
        else:
            blocklist_num = len(dataf[dataf[col_id].isin(blocklist["ID"])])
            if blocklist_num == 0:
                print('There is no matching sequences to block, check that the format of the blocklist and fasta file have the same format')
            dataf = dataf[~dataf[col_id].isin(blocklist["ID"])]
        return dataf,blocklist_num


    if args.blocklist:
        blocklist = pd.read_csv(args.blocklist, names=["ID"])
        sequences, block_num = filter_function(args.blocklist,sequences,args.gis_id,args.allowlist)
    for ch in args.delete:
        sequences["Sequence"] = sequences["Sequence"].str.replace(ch, "-", regex=False, case=False)

    with open(args.output, "w") as ofh:
        for index, row in sequences.iterrows():
            ofh.write(">{}\n{}\n".format(row["ID"], row["Sequence"]))
    print("{} sequences allowed/blacklisted from {} in the original {} file".format(block_num,len(blocklist['ID']),args.blocklist))


if __name__ == "__main__":
    main()
