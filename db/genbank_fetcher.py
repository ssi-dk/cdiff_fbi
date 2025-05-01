import argparse
from Bio import Entrez, SeqIO

def fetch_and_print_genbank_features(email: str, accession: str, locus_filter=None):
    """Fetch a GenBank record and print CDS feature info. Match by gene, product, or locus tags."""
    Entrez.email = "Your.Name.Here@example.org" #email
    print(f"Fetching GenBank record for {accession} from NCBI...")

    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")
    except Exception as e:
        print(f"Failed to fetch or parse {accession}: {e}")
        return

    print(f"\nRecord ID: {record.id}")
    print(f"Description: {record.description}")
    print(f"Sequence Length: {len(record.seq)} bp")

    print("\nCoding sequence (CDS) coordinates:")
    found = False
    matched_loci = set()

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene_name = feature.qualifiers.get("gene", [""])[0] or "unknown"
        product = feature.qualifiers.get("product", [""])[0]
        locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
        old_locus_tag = feature.qualifiers.get("old_locus_tag", [""])[0]
        gene_synonyms = feature.qualifiers.get("gene_synonym", [""])
        note = feature.qualifiers.get("note", [""])[0]

        matched = False
        comments = []
        match_sources = {}

        if locus_filter:
            for target in locus_filter:
                matched_fields = []

                if gene_name == target:
                    matched_fields.append(f"Matched by gene ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["gene"]
                    matched_loci.add(target)

                if target.lower() in product.lower():
                    matched_fields.append(f"Matched by product (partial) ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["product"]
                    matched_loci.add(target)

                if locus_tag == target:
                    matched_fields.append(f"Matched by locus_tag ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["locus_tag"]
                    matched_loci.add(target)

                if old_locus_tag == target:
                    matched_fields.append(f"Matched by old_locus_tag ({target})")
                    match_sources[target] = match_sources.get(target, []) + ["old_locus_tag"]
                    matched_loci.add(target)

                if matched_fields:
                    matched = True
                    comments.extend(matched_fields)

            # Print warning if multiple fields match the same target
            for target, fields in match_sources.items():
                if len(fields) > 1:
                    print(f"Warning: Multiple matches for CDS '{gene_name}' on locus '{target}': {', '.join(fields)}")

        else:
            matched = True  # No filter means include all features

        if not matched:
            continue

        start = int(feature.location.start) + 1  # Convert to 1-based
        end = int(feature.location.end)
        length = end - start + 1
        strand = "+" if feature.location.strand == 1 else "-"

        print(f" - {gene_name or 'unknown'}:")
        print(f"     coordinates: {start}..{end}")
        print(f"     length: {length}")
        print(f"     strand: ({strand})")
        print(f"     Locus tag: {locus_tag}")
        print(f"     Old locus tag: {old_locus_tag}")
        print(f"     Gene synonym(s): {', '.join(gene_synonyms)}")
        print(f"     Product: {product}")
        print(f"     Note: {note}")
        if comments:
            print(f"     Comments: {', '.join(comments)}")
        found = True

    # Print unmatched search terms
    if locus_filter:
        unmatched = set(locus_filter) - matched_loci
        for missing in sorted(unmatched):
            print(f"\n - Gene name '{missing}' not found in any CDS feature.")

def main():
    parser = argparse.ArgumentParser(description="Fetch and print CDS features from a GenBank record.")
    parser.add_argument("-e", "--email", required=True, help="Your email address (required by NCBI)")
    parser.add_argument("-a", "--accession", required=True, nargs="+", help="One or more GenBank accession numbers")
    parser.add_argument("-l", "--locus", nargs="*", help="Filter by gene/locus names (e.g., tcdA CD0660 cdtA)")

    args = parser.parse_args()

    for acc in args.accession:
        fetch_and_print_genbank_features(args.email, acc, args.locus)

if __name__ == "__main__":
    main()

# python genbank_fetcher.py -e rahenriksen@gmail.com -a AM180355.1 --locus tcdA tcdC cdtA tcdB CD630_06600 CD0660 cdtB cdtS