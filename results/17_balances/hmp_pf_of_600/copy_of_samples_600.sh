#!/bin/bash

# copy all hmp samples that are either oral or fecal into a separate dir. for simplicity. 
# no need to filter files in the input later...

IN="meta_9194_regions_of_clean_600.csv"
mkdir -p samples_jplace_of_600

for line in `cat ${IN}` ; do
	name=`echo $line | cut -d"," -f 1`
	cp samples_jplace/${name}.jplace samples_jplace_of_600
done
