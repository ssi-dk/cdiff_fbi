#!/usr/bin/env python3
## summarize_gene_coverage_percentages.py

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize gene coverage (positive/negative + percentage details)")
    parser.add_argument("-c", "--covfile", required=True, help="Coverage file")
    parser.add_argument("-b", "--intervalsbed", required=True, help="BED file with gene intervals")
    parser.add_argument("-o", "--output", required=True, help="Output TSV filename")
    return parser.parse_args()

# Gene-specific coverage thresholds
gene_thresholds = {
    "tcdA": 0.9,
    "tcdB": 0.9,
    "tcdC": 0.9,
    "cdtA": 0.9,
    "cdtB": 0.9,
    "other": 0.9
}

# Desired gene column order
gene_order = ["tcdA", "tcdB", "tcdC", "cdtA", "cdtB"]

def load_coverage(covfile):
    """Load coverage values as a dict: {position: total_depth}"""
    coverage = {}
    with open(covfile, "r") as f:
        for line in f:
            if line.startswith("Locus"):
                continue
            locus, total_depth, *_ = line.strip().split(",")
            pos = int(locus.split(":")[1])
            coverage[pos] = int(total_depth)
    return coverage

def evaluate_gene(coverage, start, end, threshold):
    """Return 'positive'/'negative' and float percentage value"""
    covered = sum(1 for pos in range(start + 1, end + 1) if coverage.get(pos, 0) > 0)
    length = end - start
    frac = covered / length if length else 0
    status = "positive" if frac >= threshold else "negative"
    percent = frac * 100
    return status, percent

def main():
    args = parse_args()
    coverage = load_coverage(args.covfile)

    gene_flags = {}      # "positive" / "negative"
    gene_percents = {}   # Float percentages

    with open(args.intervalsbed, "r") as bedfile:
        for line in bedfile:
            ref, start, end, gene = line.strip().split()
            if gene not in gene_order:
                continue  # skip unexpected genes
            start, end = int(start), int(end)
            threshold = gene_thresholds[gene]
            status, percent = evaluate_gene(coverage, start, end, threshold)
            gene_flags[gene] = status
            gene_percents[gene] = percent
    # Format toxin_details
    toxin_details = ";".join([
        f"{gene}_{gene_percents[gene]:.2f}" for gene in gene_order if gene in gene_percents
    ])

    # Create output dataframe
    df = pd.DataFrame([{**gene_flags, "toxin_details": toxin_details}])
    final_columns = gene_order + ["toxin_details"]
    df = df[[col for col in final_columns if col in df.columns]]

    # Save TSV
    df.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
# python Toxin.py -b intervals.bed -c ../../test/test_cdiff_single___2405W4378/sp_cdiff_fbi/test_cdiff_single___2405W4378.coverage -o ../toxin.tsv