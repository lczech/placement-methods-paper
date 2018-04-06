#!/bin/bash
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -q sandy.q                        # Select the queue
#$ -pe mvapich12 12                  # Set the parallel environment
#$ -l h_rt=24:00:00                  # Request the time for the job

source /etc/profile.d/modules.sh
module unload gcc
module load gcc/4.9.2

SAMPLE_NAME=REPLACE_SAMPLE_NAME

BASEDIR=bacterial_vaginosis
DATADIR=${BASEDIR}/00_data

ALI=${DATADIR}/Ref.phylip
SAMPLES=${DATADIR}/queries/${SAMPLE_NAME}.fasta
TREE=${BASEDIR}/01_treesearch/best_tree.newick

PAPARA=${BASEDIR}/software/papara_nt/papara

${PAPARA} -t ${TREE} -s ${ALI} -q ${SAMPLES} -j 12 -r -n best_queries_${SAMPLE_NAME}

