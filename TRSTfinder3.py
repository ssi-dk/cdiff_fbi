#!/usr/bin/env python3
## TRSTfinder3.py
## Cdiff related analysis
## Finds repeat sequences in contigs and assigns them an id according to a db.


import re
import os
import argparse


parser = argparse.ArgumentParser(description="TRSTfinder3.py")
parser.add_argument("-i", "--contigs")
parser.add_argument("-db", "--trstdb")
parser.add_argument("-o", "--outfile")
args = parser.parse_args()


## Functions
def parse_repeat_sequences(db, file):
    """Parse repeat sequences placed in a db folder
    """
    repeat_seqs_file = os.path.join(db, file)
    fragments = dict()
    seq = ""
    name = ""
    with open(repeat_seqs_file) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    fragments[name] = seq
                name = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
        if seq:
            fragments[name] = seq
    return fragments


def parse_types(db, file, fragments):
    """Parse fragment types
    """
    types_file = os.path.join(db, file)
    types = dict()
    with open(types_file) as fh:
        for line in fh:
            try:
                key, pattern = line.split(",\t")
            except ValueError:
                continue
            seq = ""
            try:
                for entry in pattern.strip().split("-"):
                    if seq:
                        seq += fragments[entry]
                    else:
                        seq = fragments[entry]
                types[key] = seq
            except KeyError:
                print("Can't find entry:", entry)
                print(list(fragments.keys()))
                exit(1)
    return types


def revtrans(seq):
    """Translation table
    """
    trans_table = str.maketrans("ATGCatgc", "TACGtacg")
    return seq.translate(trans_table)[::-1]


## Main
# Tandem repeat loci TR6 and TR10
TR6fragments = parse_repeat_sequences(args.trstdb, "TR6_repeat_sequences.ashx")
TR10fragments = parse_repeat_sequences(args.trstdb, "TR10_repeat_sequences.ashx")
TR6 = parse_types(args.trstdb, "TR6_types.ashx", TR6fragments)
TR10 = parse_types(args.trstdb, "TR10_types.ashx", TR10fragments)


# Search FSA file for hits
with open(args.contigs, "r") as fsafile:
    header = fsafile.readline()
    sequence = "".join(fsafile.readlines()).replace("\n", "")
    rTR6 = []
    rTR10 = []
    for k, v in TR6.items():
        if re.search(v, sequence, re.IGNORECASE) or re.search(revtrans(v), sequence, re.IGNORECASE):
            rTR6.append(k)
    for k, v in TR10.items():
        if re.search(v, sequence, re.IGNORECASE) or re.search(revtrans(v), sequence, re.IGNORECASE):
            rTR10.append(k)


# Print hits to outfile
with open(args.outfile, 'w') as outfile:
    print("TRST results", file=outfile)
    print(rTR6, file=outfile)
    print(rTR10, file=outfile)
    with open(os.path.join(args.trstdb, "TRST_types.ashx")) as fh:
        hitflag = "off"
        for line in fh:
            stcomb, stA, stB = line.strip().split("\t")
            if stA in rTR6 and stB in rTR10:
                print("{}\t{}\t{}".format(stcomb, stA, stB), file=outfile)
                hitflag = "on"

    if hitflag == "off":
        print("trunknown", file=outfile)
