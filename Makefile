# Makefile

.PHONY: test update_db type process summarize

clean:	
	rm -rf test/ERR142064/

# Combined test target to run all steps
# test: serotyping processing summarizing
test: serotyping md5serotypecheck processing md5processcheck

# Step 1: Run the typing process
serotyping:
	@BIN_PATH=$$(dirname $$(which python)) && \
	echo "Conda environment bin path: $$BIN_PATH" && \
	bash cdifftyping.sh -i ERR142064 -R1 test/ERR142064_1_subset.fastq.gz -R2 test/ERR142064_2_subset.fastq.gz -c test/ERR142064_subset.fasta -o test -db db -update no

# Step 2: Check MD5 checksum
md5serotypecheck:
	cd test && \
	cat ERR142064/sp_cdiff_fbi/ERR142064.indel.vcf|grep -v '#'> ERR142064/sp_cdiff_fbi/ERR142064.indel.noheader.vcf && \
	cat ERR142064/sp_cdiff_fbi/ERR142064.snp.vcf|grep -v '#'> ERR142064/sp_cdiff_fbi/ERR142064.snp.noheader.vcf && \
	cat ERR142064/sp_cdiff_fbi/ERR142064.snp_indel.vcf|grep -v '#'>ERR142064/sp_cdiff_fbi/ERR142064.snp_indel.noheader.vcf && \
	md5sum -c Serotypesums.md5

# Step 3: Run the post-processing
processing:
	bash postcdifftyping.sh -i ERR142064 -d test -stbit "STNA;NA:NA"

# Step 4: Check MD5 checksum - but remove variable date to accurately represent MD5sum
md5processcheck:
	cd test && \
	md5sum -c Processsums.md5