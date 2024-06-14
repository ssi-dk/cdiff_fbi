#!/usr/bin/env python3
## parse_report_file.py
## Cdiff related QC analysis
## Parses a report file to generate a .csv file and a .json file.g


import argparse
import os
import re
import json


parser = argparse.ArgumentParser(description="parse_report_file.py")
parser.add_argument("-r", "--report_file", help="report_file.txt")
parser.add_argument("-w", "--wgsnumber", help="wgsnumber", default="NA")
parser.add_argument("-s", "--stbit", help="ST string", default="ST;NA:NA")
parser.add_argument("-o", "--output_dir", default="output")
args = parser.parse_args()


## Functions
def print_header_to_output(csv_outfile):
	"""Print header to a file
	"""
	header = "Name;cdtA/B;tcdA;tcdB;tcdClength;117del;A117T;TRST;TR6;TR10;ST;STalleles;WGS;tcdA:tcdB:tcdC:cdtA:cdtB"
	csv_outfile=open(csv_outfile, "w")
	print(header, file=csv_outfile)
	return header


def indellength(info):
	"""Get indel length
	"""
	try:
		# Split the info string into elements
		elements = info.split("),")
		results = []
		# Regex pattern to match the nucleotide sequences
		pattern = re.compile(r"\{?\s?\d+:\s\(\"([ATCG]+)\"\,\s\"([ATCG]+)")
		for case in elements:
			match = pattern.search(case)
			if match:
				seq1, seq2 = match.groups()
				length_diff = len(seq1) - len(seq2)
				results.append(str(length_diff))
		if not results:
			results = ["0"]
		return "_".join(results)
	except Exception as e:
		# Log the exception for debugging purposes
		print(f"An error occurred: {e}")
		return "0"


def get_sampleid(line, csv_data):
	"""Line
	Strain: cdiff1
	"""
	sampleid = "-"
	if line.startswith("Strain"):
		sampleid = line.split(":")[1].strip()
		csv_data["Name"] = sampleid
	return csv_data

	
def process_tcdAB_genes(line, csv_data):
	"""Line
	Found in tcdA:Gene is not present;64/8133;{}"""
	genes = ["tcdA", "tcdB"]
	for gene in genes:
		if line.startswith(f"Found in {gene}:"):
			answer = line.strip().split(":")[1].split(";")
			if answer[0] == "Gene is present":
				csv_data[f"{gene}"] = "+"
			csv_data["cov_info"][f"{gene}"] = answer[1]

	return csv_data


def process_tcdC_gene(line, csv_data):
	"""Line
	Found in tcdC:Gene is not present;0/700;{} 
	"""
	if line.startswith("Found in tcdC:"):
		csv_data["cov_info"][f"tcdC"] = line.split(":", 1)[1].split(";")[1]
		csv_data["tcdClength"] = indellength(line.split(":", 1)[1].split(";")[2])
		
		if "18499: ('CT', 'C')" in line:
			csv_data["117del"] = "+"
			
			remaining = csv_data["tcdClength"].split("_")
			first = remaining.pop(0)
			if first in {"4", "18", "54"}:
				csv_data["tcdClength"] = first
			elif remaining:
				csv_data["tcdClength"] = str(remaining)
			else:
				csv_data["tcdClength"] = "0"
		else:
			csv_data["117del"] = "-"
	return csv_data


def process_A117T(line, csv_data):
	"""Line
	gi|126697566|ref|NC_009089.1|
	"""
	if line.startswith("gi|"):
		print(line)
		answer = line.split()
		A117T = answer[4]
		if A117T == "A" and not line[6] == "LowQual":
			header = answer[8].split(":")
			index = header.index("DP")
			depth = answer[9].split(":")[index]
			if int(depth) > 20:
				csv_data["A117T"] = "+"
	return csv_data


def process_cdtAB_genes(line, csv_data):
	"""Line
	Found in ..
	"""
	genes = ["cdtA", "cdtB"]
	for gene in genes:
		if line.startswith(f"Found in {gene}:"):
			answer = line.split(":", 1)[1].split(";")
			if answer[0] == "Gene is present":
				csv_data[f"{gene}"] = "+"
			csv_data["cov_info"][f"{gene}"] = answer[1]
	return csv_data


def process_trst_results(line, csv_data, report):
	"""Line
	[]
	"""
	if line.startswith("TRST results"):
		next_line = next(report).strip()
		TR6 = TR10 = ""
		if next_line.startswith("["):
			TR6=next_line[1:-1].split(', ')
			TR6=','.join([element.strip("'") for element in TR6])
			csv_data["TR6"] = TR6
			TR10 = next(report).strip()[1:-1].strip("'")
			csv_data["TR10"] = TR10
			csv_data["TRST"] = next(report).split()[0].strip("'")
		if not TR6:
			csv_data["TR6"] = "Unknown"
		if not TR10:
			csv_data["TR10"] = "Unknown"
		if csv_data["TRST"] == "trunknown":
			csv_data["TRST"] = "Unknown"
	return csv_data


def parse_report(report_file, stbit, wgsnumber):
	"""Parse a text file
	"""
	with open(report_file, "r") as report:
		csv_data = {
			"Name": "-",
			"cdtA": "-",
			"cdtB": "-",
			"tcdA": "-",
			"tcdB": "-",
			"tcdClength": "0",
			"117del": "-",
			"A117T": "-",
			"TRST": "-",
			"TR6": "-",
			"TR10": "-",
			"ST": f"{stbit}",
			"WGS": f"{wgsnumber}",
			"cov_info": {"tcdA": "-", "tcdB": "-", "tcdC": "-", "cdtA": "-", "cdtB": "-"}
		}
		for line in report:
			csv_data = get_sampleid(line, csv_data)
			csv_data = process_trst_results(line, csv_data, report)
			csv_data = process_A117T(line, csv_data)
			csv_data = process_tcdAB_genes(line, csv_data)
			csv_data = process_tcdC_gene(line, csv_data)
			csv_data = process_cdtAB_genes(line, csv_data)
	return csv_data
			
			
def write_to_csv(csv_data, csv_outfile):
	"""Write to csv file
	"""
	with open(csv_outfile, "a") as csvout:
		if csv_data:
			fields = [
			# header: "Name;cdtA/B;tcdA;tcdB;tcdClength;117del;A117T;TRST;TR6;TR10;ST;STalleles;WGS;tcdA:tcdB:tcdC:cdtA:cdtB"
				csv_data["Name"],
				f"{csv_data["cdtA"]}/{csv_data["cdtB"]}",
				csv_data["tcdA"],
				csv_data["tcdB"],
				csv_data["tcdClength"],
				csv_data["117del"],
				csv_data["A117T"],
				csv_data["TRST"],
				csv_data["TR6"],
				csv_data["TR10"],
				csv_data["ST"],
				csv_data["WGS"],
				f"{csv_data["cov_info"]["tcdA"]}:{csv_data["cov_info"]["tcdB"]}:{csv_data["cov_info"]["tcdC"]}:{csv_data["cov_info"]["cdtA"]}:{csv_data["cov_info"]["cdtB"]}"
			]
			print(";".join(fields), file=csvout)


def write_to_json(csv_data, json_outfile):
	"""Write to json file
	"""
	with open(json_outfile, "w") as jsonout:
		json.dump(csv_data, jsonout)

## Main
csv_data = parse_report(args.report_file, args.stbit, args.wgsnumber)
sampleid = csv_data["Name"]
csv_outfile=os.path.join(args.output_dir, f"{sampleid}.csv")
header = print_header_to_output(csv_outfile)
write_to_csv(csv_data, csv_outfile)
json_outfile=os.path.join(args.output_dir, f"{sampleid}.json")
write_to_json(csv_data, json_outfile)
