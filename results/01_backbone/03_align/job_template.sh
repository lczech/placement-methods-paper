#!/bin/bash

DATASET="hmp"
REFERENCE="Bacteria"
REFSHORT="Bact"
CONSTRAINT="Constr"

module load gcc/4.9.3
module load boost/1.60
module load openmpi/gcc

BASEDIR=/path/to/data/
NUM_TASKS=28

export OMP_NUM_THREADS=28
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

BATCH=REPLACE_BATCH

echo "This is the job for batch ${BATCH}"
echo "Running on `hostname`"
echo "Start at `date`"
echo

ALI=${BASEDIR}/00_reference/${REFERENCE}/Ref.phylip
TREE=${BASEDIR}/01_trees/${REFSHORT}_${CONSTRAINT}_backbone/best_tree.newick

QUERYDIR=${BASEDIR}/02_sequences/${DATASET}/chunks
QUERIES=REPLACE_QUERIES

PAPARA=${BASEDIR}/software/papara_nt-mpi/papara_mpi

# Run Forrest, run!
echo "Running papara..."
mpiexec -n 20 ${PAPARA} -t ${TREE} -s ${ALI} -q ${QUERIES} -j ${NUM_TASKS} -r -n ${BATCH}
echo "Finished papara"

echo
echo "Compressing and deleting..."
cd ${BASEDIR}/03_align/${DATASET}_${REFSHORT}_${CONSTRAINT}_backbone/jobs/
tar -czvf ${BATCH}.tar.gz ${BATCH}
rm -v ${BATCH}
echo "Finished compressing and deleting"

echo
echo "Batch job successful."
echo "End at `date`"
