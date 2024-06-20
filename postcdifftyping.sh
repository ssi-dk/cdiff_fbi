#!/bin/bash
## postcdifftyping.sh
## Cdiff related analysis
## Parses difftyping.sh results to generate an assessment report.


# Exit immediately if a command exits with a non-zero status.
set -e


## Functions
process_inputs() {
	while [[ $# -gt 0 ]]; do    # process all args
		case "$1" in    # assigns the args to each flag
			-i)
				sampleid="$2"
				shift 2
				;;
			-d)
				indir="$2"
				shift 2
				;;
            -stbit)
				stbit="$2"
				shift 2
				;;
            *)
				echo "Invalid option $1" >&2
				echo "Usage: $0 -i <sampleid> -d <inputdir> -stbit <stbit>"
				exit 1
				;;
		esac
	done

    # Check if all required arguments are provided
	if [ -z "$sampleid" ] ||[ -z "$indir" ] || [ -z "$stbit" ] ; then
		echo "Missing arguments!"
		echo "Usage: $0 -i <sampleid> -d <indir> -stbit <stbit>"
		exit 1
	fi
}


# Call the function with the arguments
process_inputs "$@"
echo -e "\n# Running..." 
echo "$0" "$@"


## Main
# Hardcoded paths
spcdifffbidir=$indir/$sampleid/sp_cdiff_fbi  # cdifftyping.sh results
prefix="$spcdifffbidir/$sampleid"  # prefix for indexes
rundir=$(basename "$indir")
echo "Rundir: $rundir"
wgsnumber=$(echo "$rundir"| grep -oE "N_WGS_[0-9]{3}") || wgsnumber=$rundir
echo "Wgsnumber: $wgsnumber"
report=$spcdifffbidir/${sampleid}_assessment_report.txt
echo "Report: $report"


# Making assessment report
echo -e "\n# Making assessment report..."
if [ -e $report ]; then
	echo "Skipping making assessment report... report exists: $report"
else
	echo "****Assessment report****
Strain: $sampleid" > $report

# Extracting SNP of interest at position 18500 with bash
echo -e "\n# Extracting SNP of interest at position 18500..."
grep "\s18500" "$prefix.snp.vcf" >> "$report" || true

# Extracting gene presence/absence from info files with bash
echo -e "\n# Extracting gene presence/absence from info files with bash..."
deletionA=$(head -n 1 ${prefix}_tcdA.info)
deletionB=$(head -n 1 ${prefix}_tcdB.info)
deletionC=$(head -n 1 ${prefix}_tcdC.info)
deletioncdtA=$(head -n 1 ${prefix}_cdtA.info)
deletioncdtB=$(head -n 1 ${prefix}_cdtB.info)

# Writing genes presence/absence to report with bash
echo -e "\n# Writing genes presence/absence to report with bash..."
echo "Found in tcdA:$deletionA " >> $report
echo "Found in tcdB:$deletionB " >> $report
echo "Found in tcdC:$deletionC " >> $report
echo "Found in cdtA:$deletioncdtA " >> $report
echo "Found in cdtB:$deletioncdtB " >> $report

# Extracting trst results with bash
echo -e "\n# Extracting gene presence/absence from info files with bash..."
cat "$prefix"_TRST.fasta"" >> $report
fi


# Generating qcfile with python
echo -e "\n# Generating qcfile with python..."
cmd="python3 parse_report_file.py -r $report -w $wgsnumber -s '$stbit' -o $indir/$sampleid"
if [ -e "$indir/$sampleid/$sampleid.csv" ] && [ -e "$indir/$sampleid/$sampleid.json" ]; then
	echo "Skipping generating qcfile... qcfiles exist: $indir/$sampleid/$sampleid.csv, $indir/$sampleid/$sampleid.json"
else
	echo $cmd
	eval $cmd
fi
