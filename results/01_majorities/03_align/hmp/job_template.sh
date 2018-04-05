#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 320
#$ -binding linear:12
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00


module load gcc/4.9
module load boost/1.56_gcc
module load mpi.intel/5.1_gcc

BASEDIR=/path/to/data/

NUM_TASKS=28
export OMP_NUM_THREADS=28
export MP_SINGLE_THREAD=no
export MP_TAS

BATCH=REPLACE_BATCH

echo "This is the job for batch ${BATCH}"
echo "Running on `hostname`"
echo "Start at `date`"
echo

ALI=${BASEDIR}/01_majorities/00_reference/Bacteria_sequences.fasta.reduced
TREE=${BASEDIR}/01_majorities/01_tree/best_tree.newick

QUERYDIR=${BASEDIR}/01_majorities/02_sequences/hmp/chunks
QUERIES=REPLACE_QUERIES

PAPARA=${BASEDIR}/software/papara_nt-mpi/papara_mpi

# Run Forrest, run!
echo "Running papara..."
mpiexec -n 20 ${PAPARA} -t ${TREE} -s ${ALI} -q ${QUERIES} -j ${NUM_TASKS} -r -n ${BATCH}
echo "Finished papara"

echo
echo "Compressing and deleting..."
cd ${BASEDIR}/01_majorities/03_align/hmp/jobs/
tar -czvf ${BATCH}.tar.gz ${BATCH}
rm -v ${BATCH}
echo "Finished compressing and deleting"

echo
echo "Batch job successful."
echo "End at `date`"
