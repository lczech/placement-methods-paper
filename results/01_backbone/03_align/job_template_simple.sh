#!/bin/bash

REFERENCE="Bacteria"
REFSHORT="Bact"
CONSTRAINT="Constr"

module load gcc/4.9.3
module load boost/1.60
module load openmpi/gcc

BASEDIR=/path/to/data/
NUM_TASKS=12

CHUNK=0

echo "This is the job for chunk ${CHUNK}"
echo "Start at `date`"
echo

ALI=${BASEDIR}/00_reference/${REFERENCE}/Ref.phylip
TREE=${BASEDIR}/01_trees/${REFSHORT}_${CONSTRAINT}_backbone/best_tree.newick
QUERY=${BASEDIR}/02_sequences/bv/chunks/chunk_${CHUNK}.fasta

PAPARA=${BASEDIR}/software/papara_nt-master/papara

# Run Forrest, run!
${PAPARA} -t ${TREE} -s ${ALI} -q ${QUERY} -j ${NUM_TASKS} -r -n ${CHUNK}

echo
echo "Chunk job successful."
echo "End at `date`"

