#!/bin/bash

SAMPLES="/03_bv/03_epa/orig_queries_jplace_clean/"
OUT="03_bv/08_squash_kmeans"

PROG="genesis/bin/apps/jplace_squash_kmeans"

mkdir -p ${OUT}
echo `date`

${PROG} 4 3 ${SAMPLES} ${OUT}

echo `date`
