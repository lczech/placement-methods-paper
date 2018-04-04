#!/bin/bash

for file in `cat process_scripts/run_accessions`; do
	s1=`stat fastq_gz/${file}_1.fastq.gz --printf="%s"`
	s2=`stat fastq_gz/${file}_2.fastq.gz --printf="%s"`

	sum=$(( ${s1} + ${s2} ))

	sf=`stat fasta/${file}.fasta --printf="%s"`

	echo "${file}	${sum}	${sf}"
done

