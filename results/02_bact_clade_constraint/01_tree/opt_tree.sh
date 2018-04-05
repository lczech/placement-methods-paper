#!/bin/bash

NUM_TASKS=16

BASEDIR=/path/to/here/

TREE=${BASEDIR}/01_tree/best_tree.newick
ALI=${BASEDIR}/00_reference/tax_cons_border.fasta.reduced

RAXML=${BASEDIR}/../software/standard-RAxML-master/raxmlHPC-PTHREADS-AVX

${RAXML} -f e -t ${TREE} -m GTRGAMMA -s ${ALI} -n opt_tree -T ${NUM_TASKS}
