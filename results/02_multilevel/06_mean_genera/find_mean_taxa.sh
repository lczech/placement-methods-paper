#!/bin/bash

rm -f taxa_*

for genus in `cat genera` ; do
	echo "Genus ${genus}"
	
	grep -ri "_${genus}" subtree_taxa/* >> taxa_details
	
done

wc -l taxa_details
cat taxa_details | sed 's/.*:\(.*\)\t.*/\1/g' | sort -u > taxa_filtered
wc -l taxa_filtered
wc -l taxa_filtered >> taxa_counts.txt
