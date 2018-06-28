#!/bin/bash

# Package
BACKBONE="Gen"
TAXCONSTR="Unconstr"
ALIGNMENT="General"

# Input & Output Paths
BASEDIR=/path/to/here

TREE=${BASEDIR}/11_single_seqs/01_trees/${BACKBONE}_${TAXCONSTR}/best_tree.newick
REF_ALI=${BASEDIR}/11_single_seqs/00_reference/${ALIGNMENT}_sequences.cleaned.fasta
QRY_ALI="${BASEDIR}/../data/silva/600k_taxa/600k_taxa.fasta.bin"

WORKDIR="${BASEDIR}/11_single_seqs/04_epa/${BACKBONE}_${TAXCONSTR}/"
EPA="${BASEDIR}/software/epa-ng/bin/epa-ng"

# OMP Settings
NUM_TASKS=24
export OMP_NUM_THREADS=${NUM_TASKS}
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

echo "Start at `date`"
echo "Running on `hostname`"
echo "Using ${NUM_TASKS} threads"
echo
echo "Calling: ${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --outdir ${WORKDIR} --threads ${NUM_TASKS} --opt-ref-tree --dyn-heur 0.98"
echo

# Run Forrest, run!
# cd ${BASEDIR}/software/epa/bin
${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --outdir ${WORKDIR} --threads ${NUM_TASKS} --opt-ref-tree --dyn-heur 0.98

echo
echo "End at `date`"
