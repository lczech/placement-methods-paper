#!/bin/bash

BASE="viz_correlation_trees"

SAMPLES="${BASE}/06_samples/tara_Gen_Unconstr_backbone/samples_bplace"
META="${BASE}/00_datasets/tara/meta/data_b.csv"
OUT="${BASE}/09_viz_correlation_trees/tara_Gen_Unconstr_backbone/"
#SEQS="/00_datasets/bv/meta/sequence_names.csv"

PROG="${BASE}/genesis/bin/apps/correlation_trees"

mkdir -p ${OUT}

${PROG} ${SAMPLES} ${META} ${OUT} 
#${SEQS}
