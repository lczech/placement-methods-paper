#!/bin/bash

DATASET="bv"
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

WORKDIR="${BASEDIR}/04_epa/${PACKAGE}"
EPA="EPA-ng"

echo "Workdir is ${WORKDIR}"
echo "Running on `hostname`"
echo "Start at `date`"
echo

# Input files
TREE=${BASEDIR}/01_trees/${BACKBONE}_${TAXCONSTR}_backbone/best_tree.newick
ALI=${BASEDIR}/03_align/${PACKAGE}/papara_alignment.0

# Convert to fasta
if [ ! -f ${ALI}.fasta ]; then
    echo "Converting alignments..."
    ${BASEDIR}/software/genesis/bin/phylip_fasta_conv ${ALI}
    echo "Finished converting"
    echo
fi

# Run Forrest, run!
echo "Running EPA"
cd ${BASEDIR}/software/epa/bin

hostfile=${WORKDIR}/hostfile
awk '{printf("%s\n",$1);}' $PE_HOSTFILE > ${hostfile}
mpiexec --hostfile ${hostfile} -n 1 ${EPA} -g 0.98 -t ${TREE} -s ${ALI}.fasta -O -w ${WORKDIR}

echo "Finished EPA"

echo
echo "Batch job finished."
echo "End at `date`"
