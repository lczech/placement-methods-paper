#!/bin/bash
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -q sandy.q                        # Select the queue
#$ -pe mvapich12 12                  # Set the parallel environment
#$ -l h_rt=24:00:00                  # Request the time for the job

SEED=REPLACE_SEED

BASEDIR=bacterial_vaginosis
DATADIR=${BASEDIR}/00_data
ALI=${DATADIR}/Ref.fasta

RAXML=${BASEDIR}/software/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3

${RAXML} -f o -p ${SEED} -m GTRGAMMA -s ${ALI} -n s${SEED} -T 12
