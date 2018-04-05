#!/bin/bash

# Package
BACKBONE="Bact"
TAXCONSTR="Unconstr"

# Input & Output Paths
BASEDIR=/path/to/here/
NUM_TASKS=18

TREE=${BASEDIR}/02_bact_clade_constraint/01_tree/RAxML_result.opt_tree
PARAMS=${BASEDIR}/02_bact_clade_constraint/01_tree/RAxML_info.opt_tree
REF_ALI=${BASEDIR}/02_bact_clade_constraint/00_reference/tax_cons_border_no_silva_gaps.fasta

QRY_ALI="${BASEDIR}/data/silva/600k_taxa/600k_taxa_no_silva_gaps.fasta"

WORKDIR="${BASEDIR}/02_bact_clade_constraint/02_epa/"
EPA="${BASEDIR}/software/epa-ng/bin/epa-ng"

# OMP Settings
export OMP_NUM_THREADS=${NUM_TASKS}
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

echo "Start at `date`"
echo "Running on `hostname`"
echo "Using ${NUM_TASKS} threads"
echo

# Run Forrest, run!
# cd ${BASEDIR}/software/epa/bin
${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --outdir ${WORKDIR} --threads ${NUM_TASKS} --chunk-size 10000 --opt-ref-tree

echo
echo "End at `date`"


