#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def parse_fasta(filepath):
    """Parse FASTA using Biopython's SeqIO into a dictionary."""
    return {record.id: str(record.seq) for record in SeqIO.parse(filepath, "fasta")}


def parse_types(filepath, fragments):
    """Parse TR type definitions based on repeat fragment names."""
    types = {}
    with open(filepath) as f:
        for line in f:
            if ",\t" in line:
                key, pattern = line.strip().split(",\t")
                try:
                    types[key] = ''.join(fragments[p] for p in pattern.split("-"))
                except KeyError:
                    continue
    return types


def find_matches(sequence, patterns):
    """Search for matching patterns or their reverse complements."""
    hits = []
    for name, pattern in patterns.items():
        pattern_seq = Seq(pattern)
        if re.search(str(pattern_seq), sequence, re.IGNORECASE) or re.search(str(pattern_seq.reverse_complement()), sequence, re.IGNORECASE):
            hits.append(name)
    return hits


def load_trst_types(filepath):
    """Load TRST type mapping from a TSV file."""
    trst_table = []
    with open(filepath) as f:
        for line in f:
            trst, tr6, tr10 = line.strip().split("\t")
            trst_table.append((trst, tr6, tr10))
    return trst_table


def match_trst(rTR6, rTR10, trst_table):
    """Find the matching TRST type from TR6/TR10 combination."""
    for trst, tr6, tr10 in trst_table:
        if tr6 in rTR6 and tr10 in rTR10:
            return trst
    return "Unknown"


def run_trst_typing_on_fasta(fasta_path, db_dir):
    """Run TRST typing on a FASTA file. Returns a pandas DataFrame."""
    print(f"Loading database from: {db_dir}")

    TR6_frags = parse_fasta(os.path.join(db_dir, "TR6_repeat_sequences.fa"))
    TR10_frags = parse_fasta(os.path.join(db_dir, "TR10_repeat_sequences.fa"))

    TR6_types = parse_types(os.path.join(db_dir, "TR6_types.txt"), TR6_frags)
    TR10_types = parse_types(os.path.join(db_dir, "TR10_types.txt"), TR10_frags)

    trst_table = load_trst_types(os.path.join(db_dir, "TRST_types.txt"))
    print(f"Done loading fasta files and TRST types")

    all_TR6_hits = set()
    all_TR10_hits = set()
    contig_count = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        contig_count += 1
        seq = str(record.seq)

        rTR6 = find_matches(seq, TR6_types)
        rTR10 = find_matches(seq, TR10_types)

        all_TR6_hits.update(rTR6)
        all_TR10_hits.update(rTR10)

        if contig_count % 1000 == 0:
            print(f"Processed {contig_count} contigs...")

    print(f"Completed processing {contig_count} contigs.")

    for trst, tr6, tr10 in trst_table:
        if tr6 in all_TR6_hits and tr10 in all_TR10_hits:
            return {"TRST": trst, "TR6": tr6, "TR10": tr10}

    return {"TRST": "Unknown", "TR6": "Unknown", "TR10": "Unknown"}


def main():
    parser = argparse.ArgumentParser(description="TRST Typer - Identify TR6/TR10 patterns and assign TRST types.")
    parser.add_argument("-i", "--input_fasta", required=True, help="Input contigs FASTA file")
    parser.add_argument("-d", "--db_dir", required=True, help="Directory containing TRST DB files")
    parser.add_argument("-o", "--output", required=True, help="Output .tsv file to save results")
    args = parser.parse_args()

    print(f"Starting TRST typing on file: {args.input_fasta}")
    result = run_trst_typing_on_fasta(args.input_fasta, args.db_dir)

    with open(args.output, "w") as out:
        out.write("TRST\tTR6\tTR10\n")
        out.write(f"{result['TRST']}\t{result['TR6']}\t{result['TR10']}\n")

    print(f"TRST results saved to: {args.output}")

if __name__ == "__main__":
    main()

#python TRST/TRST_typer.py -i ../test/test_cdiff_single___2405W4378/test_cdiff_single___2405W4378.fasta -d TRST/ -o sample_TRST_typing.tsv