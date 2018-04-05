#!/bin/bash



BASEDIR="11_consensus_seqs/00_reference"

cd ${BASEDIR}
for dir in `find * -maxdepth 0 -type d`; do

	echo "================================" >> reduce_silva_gaps.log
	echo "Dir: ${dir}" >> reduce_silva_gaps.log
	cd ${dir}
	
	genesis/bin/apps/silva_gap_user Archaea_sequences.fasta    >> ../reduce_silva_gaps.log
	genesis/bin/apps/silva_gap_user Bacteria_sequences.fasta   >> ../reduce_silva_gaps.log
	genesis/bin/apps/silva_gap_user Eukaryota_sequences.fasta  >> ../reduce_silva_gaps.log
	genesis/bin/apps/silva_gap_user General_sequences.fasta    >> ../reduce_silva_gaps.log
	
	cd ..
	
done
