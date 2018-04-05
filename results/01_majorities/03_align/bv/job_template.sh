#!/bin/bash
#
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -q sandy.q                        # Select the queue
#$ -pe mvapich12 12                  # Set the parallel environment
#$ -l h_rt=24:00:00                  # Request the time for the job

module load gcc/4.9.3
module load boost/1.60
module load openmpi/gcc

BASEDIR=/path/to/data/
NUM_TASKS=12

CHUNK=0

echo "This is the job for chunk ${CHUNK}"
echo "Start at `date`"
echo

ALI=${BASEDIR}/01_majorities/00_reference/Bacteria_sequences.fasta.reduced
TREE=${BASEDIR}/01_majorities/01_tree/best_tree.newick
QUERY=${BASEDIR}/01_majorities/02_sequences/bv/chunks/chunk_${CHUNK}.fasta

PAPARA=${BASEDIR}/software/papara_nt-master/papara

# Run Forrest, run!
${PAPARA} -t ${TREE} -s ${ALI} -q ${QUERY} -j ${NUM_TASKS} -r -n ${CHUNK}

echo
echo "Chunk job successful."
echo "End at `date`"

