#!/bin/bash

SEED=REPLACE_SEED
NUM_TASKS=28

BASEDIR=/path/to/data

DATADIR=${BASEDIR}/00_reference
ALI=${DATADIR}/Bacteria/Ref.phylip

RAXML=${BASEDIR}/software/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3

${RAXML} -f o -p ${SEED} -m GTRGAMMAX -s ${ALI} -n s${SEED} -T ${NUM_TASKS}
