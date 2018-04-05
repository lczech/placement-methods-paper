#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 16
#$ -binding linear:16
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

source /etc/profile.d/modules.sh

# Package
BACKBONE="Arch"
TAXCONSTR="Constr"

# Input & Output Paths
BASEDIR=/path/to/results/

TREE=${BASEDIR}/01_trees/${BACKBONE}_${TAXCONSTR}_backbone/best_tree.newick
REF_ALI=${BASEDIR}/00_reference/${BACKBONE}/tax_cons_border.fasta
QRY_ALI="${BASEDIR}/../data/silva/600k_taxa.fasta.bin"

WORKDIR="${BASEDIR}/11_tree_eval/${BACKBONE}_${TAXCONSTR}_backbone/"
EPA="${BASEDIR}/software/epa-ng/bin/epa-ng"

# Load Modules
module unload gcc
module load gcc/4.9.3
module load binutils/2.25
module load gnutools/1.0

# OMP Settings
NUM_TASKS=16
export OMP_NUM_THREADS=${NUM_TASKS}
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

echo "Start at `date`"
echo "Running on `hostname`"
echo
echo "Calling: ${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --out-dir ${WORKDIR} --threads ${NUM_TASKS} --opt-ref-tree --dyn-heur 0.98"
echo

# Run Forrest, run!
# cd ${BASEDIR}/software/epa/bin
${EPA} --ref-msa ${REF_ALI} --tree ${TREE} --query ${QRY_ALI} --outdir ${WORKDIR} --threads ${NUM_TASKS} --opt-ref-tree --dyn-heur 0.98

echo
echo "End at `date`"
