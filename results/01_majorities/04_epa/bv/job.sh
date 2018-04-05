#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich12 12
#$ -binding linear:12
#$ -q sandy.q
#$ -cwd
#$ -j y 
#$ -l h_rt=24:00:00

module load gcc/4.9.3
module load boost/1.60
module load openmpi/gcc

BASEDIR=/path/to/data/
NUM_TASKS=12

WORKDIR="${BASEDIR}/01_majorities/04_epa/bv"
EPA=${BASEDIR}/software/epa/bin/epa

echo "Workdir is ${WORKDIR}"
echo "Running on `hostname`"
echo "Start at `date`"
echo

# Input files
TREE=${BASEDIR}/01_majorities/01_tree/best_tree.newick
ALI=${BASEDIR}/01_majorities/03_align/bv/papara_alignment.0

# Convert to fasta
if [ ! -f ${ALI}.fasta ]; then
    echo "Converting alignments..."
    ${BASEDIR}/software/genesis/bin/phylip_fasta_conv ${ALI}
    echo "Finished converting"
    echo
fi


# Run Forrest, run!
echo "Running EPA"

hostfile=${WORKDIR}/hostfile
awk '{printf("%s\n",$1);}' $PE_HOSTFILE > ${hostfile}
mpiexec --hostfile ${hostfile} -n 1 ${EPA} -g 0.99 -t ${TREE} -s ${ALI}.fasta -O -w ${WORKDIR}

echo "Finished EPA"


echo
echo "Batch job finished."
echo "End at `date`"
