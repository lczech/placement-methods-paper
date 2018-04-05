#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 16
#$ -binding linear:16
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

SEED=REPLACE_SEED

BASEDIR=/path/to/here
NUM_TASKS=28


ALI=${BASEDIR}/02_multilevel/00_reference/REPLACE_METHOD/Bacteria_REPLACE_REFERENCE.tips.reduced.fasta

RAXML=${BASEDIR}/software/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3

${RAXML} -f o -p ${SEED} -m GTRGAMMAX -s ${ALI} -n s${SEED} -T ${NUM_TASKS}
