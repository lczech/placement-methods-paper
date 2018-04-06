#!/bin/bash

SAMPLES="/03_bv/03_epa/orig_queries_jplace_clean"
META="data/bv/meta/meta_simple.csv"
SEQS="/00_datasets/bv/meta/sequence_names.csv"
OUT="/01_backbone/09_viz_correlation_trees/bv_reftree/"

PROG="genesis/bin/apps/correlation_trees"

mkdir -p ${OUT}

${PROG} ${SAMPLES} ${META} ${OUT} ${SEQS}
