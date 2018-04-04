#!/bin/bash

DATASET="hmp"
BACKBONE="Bact"
TAXCONSTR="Constr"
PACKAGE="${DATASET}_${BACKBONE}_${TAXCONSTR}_backbone"

BASEDIR=/path/to/data

module load gcc/4.9
module load boost/1.56_gcc
module load mpi.intel/5.1_gcc

NUM_TASKS=28
export OMP_NUM_THREADS=28
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

BATCH=REPLACE_BATCH
WORKDIR="${BASEDIR}/04_epa/${PACKAGE}/jobs/${BATCH}"
EPA="EPA-ng"

echo "This is the job for batch ${BATCH}"
echo "Workdir is ${WORKDIR}"
echo "Running on `hostname`"
echo "Start at `date`"
echo

# unpack alignment
echo
echo "Extracting alignments..."
cd ${BASEDIR}/03_align/${PACKAGE}/jobs/
tar xvzf ${BATCH}.tar.gz
echo "Finished extracting"
echo

# Convert to fasta
echo "Converting alignments..."
for chunk in `ls ${BATCH}/papara_alignment.*`; do
    ${BASEDIR}/software/genesis/bin/phylip_fasta_conv ${BASEDIR}/03_align/${PACKAGE}/jobs/${chunk}
done
echo "Finished converting"
echo

# Chunk list
cd ${BASEDIR}/03_align/${PACKAGE}/jobs/${BATCH}
ls papara_alignment.*.fasta > chunk_list
QUERYDIR=${BASEDIR}/03_align/${PACKAGE}/jobs/${BATCH}
#QUERIES=`awk -vORS=, '{ print "${QUERYDIR}/"$1 }' chunk_list | sed 's/,$/\n/'`
QUERIES=`awk -v QUERYDIR="${QUERYDIR}" -vORS=, '{ print QUERYDIR"/"$1 }' chunk_list | sed 's/,$/\n/'`
echo "queries: ${QUERIES}"
echo

# Input files
TREE=${BASEDIR}/01_trees/${BACKBONE}_${TAXCONSTR}_backbone/best_tree.newick
# ALI=${BASEDIR}/03_align/${PACKAGE}/jobs/${BATCH}/papara_alignment.${CHUNK}.fasta

# Run Forrest, run!
echo "Running EPA"
cd ${BASEDIR}/software/epa/bin

hostfile=${WORKDIR}/hostfile
awk '{printf("%s\n",$1);}' $PE_HOSTFILE > ${hostfile}
mpiexec -n 20 ${EPA} -g 0.98 -t ${TREE} -s ${QUERIES} -O -w ${WORKDIR}

echo "Finished EPA"

echo
echo "Files in align dir ${BASEDIR}/03_align/${PACKAGE}/jobs/${BATCH}"
ls ${BASEDIR}/03_align/${PACKAGE}/jobs/${BATCH}

echo
echo "Files in workdir ${WORKDIR}"
ls ${WORKDIR}

# deleting alignment
echo
echo "Deleting alignment..."
cd ${BASEDIR}/03_align/${PACKAGE}/jobs/
rm -r ${BATCH}/
echo "Finished deleting"

echo
echo "Batch job finished."
echo "End at `date`"
