#!/bin/bash



BASEDIR="path/to/02_multilevel/02_sequences/"

cd ${BASEDIR}
for file in `ls *.fasta`; do

	echo "================================" >> reduce_silva_gaps.log
	echo "File: ${file}" >> reduce_silva_gaps.log
	
	genesis/bin/apps/silva_gap_user ${file}    >> reduce_silva_gaps.log
	
done
