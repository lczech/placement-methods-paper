#!/bin/bash

rm -f genus_*
rm -f species_*
rm -f species_counts.txt
rm -f genus_counts.txt

mkdir -p genus_details
mkdir -p species_details

rm -f genus_details/*
rm -f species_details/*

for species in `cat species` ; do
	filename="${species/\.\*/_}"
	echo "Species ${filename}"
	
	egrep -ri "_${species}" ../04_epa/cavener/*_Unconstr | grep "SEQ_" >> species_details/${filename}
	wc -l species_details/${filename}
	
	cat species_details/${filename} | sed 's/.*\(SEQ_[0-9]*_\).*/\1/g' | sort -u > species_details/${filename}_filtered
	wc -l species_details/${filename}_filtered
	wc -l species_details/${filename}_filtered >> species_counts.txt
done

for genus in `cat genera` ; do
	echo "Genus ${genus}"
	
	grep -ri "_${genus}" ../04_epa/cavener/*_Unconstr | grep "SEQ_" >> genus_details/${genus}
	wc -l genus_details/${genus}
	
	cat genus_details/${genus} | sed 's/.*\(SEQ_[0-9]*_\).*/\1/g' | sort -u > genus_details/${genus}_filtered
	wc -l genus_details/${genus}_filtered
	wc -l genus_details/${genus}_filtered >> genus_counts.txt
done


