#!/bin/bash

let total=0
let filtered=0

for i in {0..63}; do
echo "Sample ${i}"
let subtotal=$(grep ">" 2017.12.29_11.37.26_sample_${i}/reads/sample_${i}.fasta | wc -l )
let total+=$subtotal
let subfiltered=$(grep ">" 2017.12.29_11.37.26_sample_${i}/reads/sample_${i}.1* | wc -l )
let filtered+=$subfiltered
echo "  ${subfiltered} / ${subtotal}"
done
echo "Final Tally"
echo "${filtered} / ${total}"
#2017.12.29_11.37.26_sample_0/reads/sample_0.18S_rRNA
