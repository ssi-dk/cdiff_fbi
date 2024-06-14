#!/usr/bin/python3
## extract_genes_and_indels_from_vcf.py
## Cdiff related analysis
## Extracts genes and indels from a .vcf file according to intervals defined in a .bed file.


import argparse


parser = argparse.ArgumentParser(description="extract_genes_and_indels_from_vcf.py")
parser.add_argument("-i", "--indelvcf", help="sampleid.indel.vcf")
parser.add_argument("-o", "--outputname", help="outputname_gene.info")
parser.add_argument("-c", "--covfile", help="sampleid.coverage")
parser.add_argument("-b", "--intervalsbed", help="intervals.bed")
args = parser.parse_args()


## Functions
def get_indels(indelvcf, start_pos, end_pos):
	"""Get indels in a specific genomic region
	"""
	indels_dict = dict()
	with open(indelvcf, "r") as infile:
		for line in infile:
			if line.startswith("#"):
				continue
			(chrom, pos, ID, ref, alt, qual, flt, info, frmt, remainder) = line.split("\t")
			pos = int(pos)
			if pos >= start_pos and pos <=end_pos:
				#info_data=parse_vcf_info(info)
				#if int(info_data["AC"])==1 and readrem(remainder)<0.3:
				indels_dict[pos]=(ref,alt)
	return indels_dict


def get_gene_coverage(covfile, start_pos, end_pos):
	"""Locus,Total_Depth,Average_Depth_sample,Depth_for_sampleid
	gi|126697566|ref|NC_009089.1|:5551,84,84.00,84
	"""
	# list of lines within gene boundaries without coverage
	lines_without_coverage = []
	presentstat=0 # Starts at 1 because the lenght of the gene includes the start position
	with open(covfile, "r") as infile:
		for line in infile:
			if line.startswith("Locus"):
				continue
			locus = line.split(",")[0]
			pos = int(locus.split(":")[1])
			total_depth = int(line.split(",")[1])
			if int(pos)>start_pos and int(pos)<=end_pos:
				if total_depth != 0:
					presentstat+=1
				else:
					lines_without_coverage.append(line)
			if int(pos)>end_pos:
				break
	return presentstat, lines_without_coverage


## Main
with open(args.intervalsbed, "r") as infile:
	for line in infile:
		# line
		# gi|126697566|ref|NC_009089.1|   9450    17582   tcdA
		line = line.strip().split()
		ref = line[0]
		start_pos = int(line[1])
		end_pos = int(line[2])
		gene = line[3]
		indels_dict = get_indels(args.indelvcf, start_pos, end_pos)
		presentstat, lines_without_coverage = get_gene_coverage(args.covfile, start_pos, end_pos)
		length=end_pos-start_pos
		with open(f"{args.outputname}_{gene}.info", "w") as outfile:
			if presentstat/length > 0.9:
				outfile.write(";".join(["Gene is present", "/".join([str(presentstat), str(length)]), ""]))
			else:
				outfile.write(";".join(["Gene is not present", "/".join([str(presentstat), str(length)]), ""]))
			print(indels_dict, file=outfile)
			if lines_without_coverage:
					for l in lines_without_coverage:
						outfile.write(f"{l}")			  
