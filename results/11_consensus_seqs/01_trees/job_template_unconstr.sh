#!/bin/bash


SEED=REPLACE_SEED

BASEDIR=/path/to/here
NUM_TASKS=12


ALI=${BASEDIR}/11_consensus_seqs/00_reference/REPLACE_METHOD/REPLACE_REFERENCE_sequences.fasta.reduced

RAXML=${BASEDIR}/software/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3

${RAXML} -f o -p ${SEED} -m GTRGAMMAX -s ${ALI} -n s${SEED} -T ${NUM_TASKS}
