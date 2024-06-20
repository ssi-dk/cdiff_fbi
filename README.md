# cdiff_fbi
Typing of *Clostridioides difficile* isolates using NGS data (reads and contigs) based on tandem repeat loci (TR6, TR10), and toxin genes (cdtA, cdtB, tcdA, tcdB, tcdC).

## Quick start
```bash
# Type
bash cdifftyping.sh -h
# Process
bash postcdifftyping.sh -h
# Summarize
python3 qc_cdiff_summary.py -h
```

## Installation
### Source
```bash
# Clone this repo
git clone https://github.com/ssi-dk/cdiff_fbi.git
# Create an environment with the required tools with conda
conda create --name cdiff_pipeline picard gatk4 biopython ruamel.yaml kraken bwa samtools
# Activate the environment
conda activate cdiff pipeline
# Install a custom tool
git clone https://github.com/ssi-dk/serum_readfilter
cd serum_readfilter
pip install .
```

## Usage
### Example
```bash
# Download data into the test folder
mkdir -p test
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/ERR142064/ERR142064_2.fastq.gz -P test
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/ERR142064/ERR142064_1.fastq.gz -P test
touch ERR142064.fasta  # Create an empty file as a fake assembly for testing purposes
```

### Pipeline
```bash
# Build the db
bash cdifftyping.sh -db db -update yes  # WARNING: Not yet implemented for serumdb or trstdb
# Type
bash cdifftyping.sh -i ERR142064 -R1 test/ERR142064_1.fastq.gz -R2 test/ERR142064_2.fastq.gz -c test/ERR142064.fasta -qc pass -o test -db db -update no
# Process
bash postcdifftyping.sh -i ERR142064 -d test -stbit "STNA;NA:NA"
# Summarize
python3 qc_cdiff_summary.py -i test -o test
```

## Output
### .csv
```
Name;cdtA/B;tcdA;tcdB;tcdClength;117del;A117T;TRST;TR6;TR10;ST;STalleles;WGS;tcdA:tcdB:tcdC:cdtA:cdtB
ERR142064;+/+;+;+;0;+;-;Unknown;Unknown;Unknown;STNA;NA:NA;test;8119/8133:6914/7101:700/700:1389/1389:2628/2628
```
### .json
```
{"Name": "ERR142064", "cdtA": "+", "cdtB": "+", "tcdA": "+", "tcdB": "+", "tcdClength": "0", "117del": "+", "A117T": "-", "TRST": "Unknown", "TR6": "Unknown", "TR10": "Unknown", "ST": "STNA;NA:NA", "WGS": "test", "cov_info": {"tcdA": "8119/8133", "tcdB": "6914/7101", "tcdC": "700/700", "cdtA": "1389/1389", "cdtB": "2628/2628"}}
```

## Updating the db
```bash
# Build the db
bash cdifftyping.sh -db db -update yes  # WARNING: Not yet implemented for serumdb or trstdb
```
