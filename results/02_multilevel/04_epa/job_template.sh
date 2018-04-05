#!/bin/bash

source /etc/profile.d/modules.sh

# Package
METHOD="REPLACE_METHOD"
BACKBONE="REPLACE_REFERENCE"
TAXCONSTR="REPLACE_CONSTR"

# Input & Output Paths
BASEDIR=/path/to/here
NUM_TASKS=28

TREE=${BASEDIR}/02_multilevel/01_trees/${METHOD}/${BACKBONE}_${TAXCONSTR}/best_tree.newick
REF_ALI=${BASEDIR}/02_multilevel/00_reference/${METHOD}/Bacteria_${BACKBONE}.tips_no_silva_gaps.fasta
#QRY_ALI="${BASEDIR}/../data/silva/600k_taxa/600k_taxa_no_silva_gaps.fasta.bin"
QRY_ALI="${BASEDIR}/02_multilevel/02_sequences/subtree_sequences/Bacteria_${BACKBONE}_no_silva_gaps.fasta"

WORKDIR="${BASEDIR}/02_multilevel/04_epa/${METHOD}/${BACKBONE}_${TAXCONSTR}/"
# EPA="${BASEDIR}/software/epa-ng/bin/epa-ng"
EPA="${BASEDIR}/software/epa-ng-mpi/bin/epa-ng"

# Load Modules
module unload mpi.ibm/1.4
module unload gcc
module load gcc/6
module load mpi.intel/5.1_gcc
module load bison

# OMP Settings
export OMP_NUM_THREADS=${NUM_TASKS}
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

echo "Start at `date`"
echo "Running on `hostname`"
echo "Using ${NUM_NODES} nodes and ${NUM_TASKS} threads"
echo
echo "Calling: ${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --outdir ${WORKDIR} --threads ${NUM_TASKS} --opt-ref-tree --dyn-heur 0.98"
echo

# Run Forrest, run!
# cd ${BASEDIR}/software/epa/bin
${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --outdir ${WORKDIR} --threads ${NUM_TASKS} --opt-ref-tree --dyn-heur 0.98

echo
echo "End at `date`"
