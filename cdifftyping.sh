#!/bin/bash
## cdifftyping.sh
## Cdiff related analysis
## Processes sequencing data, including quality control, read alignment, variant calling (SNPs and INDELs), and gene coverage.


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
            -R1)
                read1="$2"
                shift 2
                ;;
            -R2)
                read2="$2"
                shift 2
                ;;
            -c)
                contigs="$2"
                shift 2
                ;;
            -qc)
                qcstatus="$2"
                shift 2
                ;;
            -o)
                outdir="$2"
                shift 2
                ;;
            -db)
                db="$2"
                shift 2
                ;;
            -update)
                update="$2"
                shift 2
                ;;
            *)
				echo "Unknown option: $1"
                echo "Usage: $0 [-i <sampleid>] [-R1 <read1>] [-R2 <read2>] [-c <contigs.fasta>] [-qc <qcstatus>] [-o <outdir>] -db <db> -update <update:yes/no>"
                exit 1
                ;;
        esac
    done

    # Check if the required arguments -db and -update are provided
    if [ -z "$db" ] || [ -z "$update" ]; then
        echo "Missing required arguments!"
		echo "Usage: $0 [-i <sampleid>] [-R1 <read1>] [-R2 <read2>] [-c <contigs.fasta>] [-qc <qcstatus>] [-o <outdir>] -db <db> -update <update:yes/no>"
        echo 
        exit 1
    fi
}

# Call the function with the arguments
process_inputs "$@"
echo -e "\n# Running..." 
echo "$0" "$@"


## Main
# Hardcoded paths
serumdb="$db/cdiff_serum_readfilter"
sref=$(ls $serumdb/*.fasta) # kraken.fasta
toxinsdb="$db/cdiff_toxins" 
reference=$(ls $toxinsdb/*.fasta)  # tcdRegion.fasta, NC_009089.1 (genome) 
intervals=$(ls $toxinsdb/*.bed)  # intervals.bed, toxin genes (genomic regions)
dict="$toxinsdb/$(basename "$reference" .fasta).dict"
trstdb="$db/cdiff_TRST"
echo -e "\n# References and variables..."
echo "SerumDB: $serumdb"
echo "ToxinsDB: $toxins"
echo "Reference: $reference"
head -n 1 "$reference"
echo "Intervals: $intervals"
cat "$intervals"
echo "TRSTDB: $trstdb"


# Update dbs
if [[ "$update" == "yes" ]] ; then
	echo -e "\n# Updating dbs..."
	## TODO: Find which command was used to make the serumdb, code below is a draft!
	echo -e "\n# WARNING!: Code not implemented yet!"
	echo -e "\n# Exiting.."
	exit 1
	# # cdiff_serum_readfilter
	# echo -e "\n# Updating cdiff_serum_readfilter..."
	# rm -r $serumdb/library $serumdb/taxonomy || true  # Remove old folders
	# rm $serumdb/database.* $serumdb/lca.complete $serumdb/seqid2taxid.map || true  # Remove old files
	# echo -e "\n# Building the kraken serum_readfilter db..."
	# serum_readfilter makedb kraken -db $serumdb -ref $sref

	# cdiff_toxins
	echo -e "\n# Updating cdiff_toxinsdb..."
	echo -e "\n# Building indexes with bwa and samtools..."
	bwa index $reference
	samtools faidx $reference
	echo -e "\n# Making dictionaries of the reference with picard..."
	rm $dict || true # Remove old dict
	picard CreateSequenceDictionary -R $reference -O $dict
else
	echo -e "\n# Using prebuild dbs..."
fi


# Check if any of the optional arguments were provided
if [ -n "$sampleid" ] || [ -n "$read1" ] || [ -n "$read2" ] || [ -n "$contigs" ] || [ -n "$qcstatus" ] || [ -n "$outdir" ]; then
	echo -e "\n# Arguments provided..."
		[ -n "$sampleid" ] && echo "Sample ID: $sampleid"
		[ -n "$read1" ] && echo "Read1: $read1"
		[ -n "$read2" ] && echo "Read2: $read2"
		[ -n "$contigs" ] && echo "Contigs: $contigs"
		[ -n "$qcstatus" ] && echo "QC Status: $qcstatus"
		[ -n "$outdir" ] && echo "Output Directory: $outdir"
		echo "Database: $db"
		echo "Update: $update"
else
    echo -e "\n# No further task to be done..."
    exit 1
fi


# Create outputs
rundir=$(basename $(dirname "$read1"))
echo "Rundir: $rundir"
wgsnumber=$(echo "$rundir"| grep -oE "N_WGS_[0-9]{3}") || wgsnumber=$rundir
echo "Wgsnumber: $wgsnumber"
spcdifffbidir=$outdir/$sampleid/sp_cdiff_fbi  # cdifftyping.sh results
mkdir -p $spcdifffbidir
prefix="$spcdifffbidir/$sampleid"  # prefix for indexes


# Running QC check for enough reads (> 1000K)
echo -e "\n# Running QC check for enough reads (> 1000K)..."
if [ -e $read1 ]; then
	lines=$( wc -l $read1 | awk '{print $1}' )
	if [ $lines \< 1000000 ]; then
		qcstatus="short"
		echo "Found $lines lines in $read1"
		echo "QC status: $qcstatus"
		echo "$cmd"
		eval $cmd
		exit
	else
		echo "Found $lines lines in $read1"
		echo "QC status: $qcstatus"
	fi
else
	echo "File not found: $read1"
	exit
fi


# Filtering reads with serum_readfilter
echo -e "\n# Filtering reads with serum_readfilter..."
cmd="serum_readfilter runfilter kraken -R1 ${read1} -R2 ${read2} -o $spcdifffbidir/cdifffiltered -db $serumdb"
r1="$spcdifffbidir/cdifffiltered_R1.fastq"
r2="$spcdifffbidir/cdifffiltered_R2.fastq"
filtered_reads="$r1 $r2"
if [ -e "$r1" ] && [ -e "$r2" ]; then
	echo "Skipping filtering reads... filtered reads exist: $filtered_reads"
else
	echo "$cmd"
	eval $cmd
fi


# Aligning reads to the reference using bwa
echo -e "\n# Aligning reads to the reference using bwa..."
if [ -e "$prefix.sam" ]; then
	echo "Skipping aligning reads... alignment file exists: $prefix.sam"
else
	bwa mem -t 1 $reference $filtered_reads -R "@RG\tID:1\tPL:illumina\tPU:barcode\tLB:myLibrary\tSM:$sampleid" > $prefix.sam
fi


# Translating sam to bam and sorting by position (reference) with samtools
echo -e "\n# Translating sam to bam and sorting by position (reference) with samtools..."
cmd="samtools view -b $prefix.sam | samtools sort -o $prefix.bam"
if [ -e "$prefix.bam" ]; then
	echo "Skipping translation and sorting... bam file exists: $prefix.bam "
else
	echo "$cmd"
	eval $cmd
fi


# Building indexes with picard and samtools
echo -e "\n# Building indexes with picard and samtools..."
if [ -e "$prefix.bam.bai" ]; then
	echo "Skipping translation and sorting... bam.bai file exists: $prefix.bam.bai"
else
	picard BuildBamIndex -I $prefix.bam -O $prefix.bam.bai
	samtools index $prefix.bam
fi


# Calling SNPs and INDELs with gatk
echo -e "\n# Calling SNPs and INDELs with gatk..."
if [ -e "$prefix.snp_indel.vcf" ] && [ -e "$prefix.snp.vcf" ] && [ -e "$prefix.indel.vcf" ] ; then
	echo "Skipping calling SNPs... vcf files exist: $prefix.snp_indel.vcf $prefix.snp.vcf $prefix.indel.vcf"
else
	gatk HaplotypeCaller -I $prefix.bam -R $reference -O $prefix.snp_indel.vcf -ploidy 1 --output-mode EMIT_VARIANTS_ONLY --annotate-with-num-discovered-alleles --max-alternate-alleles 6
	gatk SelectVariants -V $prefix.snp_indel.vcf -select-type SNP -O $prefix.snp.vcf
	gatk SelectVariants -V $prefix.snp_indel.vcf -select-type INDEL -O $prefix.indel.vcf
fi


# Getting gene coverage with gatk
echo -e "\n# Getting gene coverage with gatk..."
if [ -e "$prefix.coverage" ] ; then
	echo "Skipping getting gene coverage... file exists: $prefix.coverage"
else
	gatk DepthOfCoverage -R $reference -I $prefix.bam -L $intervals -O $prefix.coverage
fi


# Extracting genes and indels from vcf files with python
echo -e "\n# Extracting genes and indels from vcf files with python..."
cmd="python3 extract_genes_and_indels_from_vcf.py -i $prefix.indel.vcf -o $prefix -c $prefix.coverage -b $intervals"
echo $cmd
eval $cmd


# Finding trst with trstfinder
echo -e "\n# Finding trst with trstfinder..."
cmd="python3 TRSTfinder3.py -i $contigs -db $trstdb -o $prefix'_TRST.fasta'"
if [ -e $prefix"_TRST.fasta" ]; then
	echo "Skipping finding trst... file exists: $prefix'_TRST.fasta'"
else
	echo $cmd
	eval $cmd
fi
