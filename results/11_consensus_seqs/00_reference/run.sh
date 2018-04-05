#!/bin/bash

genesis/bin/apps/silva_consensus_seqs Archaea
genesis/bin/apps/silva_consensus_seqs Bacteria
genesis/bin/apps/silva_consensus_seqs Eukaryota
genesis/bin/apps/silva_consensus_seqs General

BASEDIR="11_consensus_seqs/00_reference"

# genesis/bin/apps/fasta_cleanup ${BASEDIR}/Archaea_sequences.fasta
# genesis/bin/apps/fasta_cleanup ${BASEDIR}/Bacteria_sequences.fasta
# genesis/bin/apps/fasta_cleanup ${BASEDIR}/Eukaryota_sequences.fasta
# genesis/bin/apps/fasta_cleanup ${BASEDIR}/General_sequences.fasta

cd ${BASEDIR}
for dir in `find * -maxdepth 0 -type d`; do

	cd ${dir}
	raxml -f c -m GTRGAMMA -s Archaea_sequences.fasta -n Archaea_reduce >> Archaea_reduce.log
	raxml -f c -m GTRGAMMA -s Bacteria_sequences.fasta -n Bacteria_reduce >> Bacteria_reduce.log
	raxml -f c -m GTRGAMMA -s Eukaryota_sequences.fasta -n Eukaryota_reduce >> Eukaryota_reduce.log
	raxml -f c -m GTRGAMMA -s General_sequences.fasta -n General_reduce >> General_reduce.log
	cd ..
	
done
