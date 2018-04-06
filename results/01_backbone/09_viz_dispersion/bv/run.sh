#!/bin/bash

SAMPLES="/03_bv/03_epa/orig_queries_jplace_clean/"
OUT="/01_backbone/09_viz_dispersion/bv/iod_no_clip/"

PROG="genesis/bin/apps/dispersion_trees"

mkdir -p ${OUT}

${PROG} ${SAMPLES} ${OUT}
