#!/bin/bash
#
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -q sandy.q                        # Select the queue
#$ -pe mvapich12 12                  # Set the parallel environment
#$ -l h_rt=24:00:00                  # Request the time for the job

BASEDIR=/path/to/results
NUM_TASKS=28

tree=Arch

iqtree=${BASEDIR}/software/iqtree-omp-1.5.6-Linux/bin/iqtree-omp
ALI=${BASEDIR}/01_backbone/13_significance_tests/${tree}/tax_cons_border.phylip
TREES=${BASEDIR}/01_backbone/13_significance_tests/${tree}/${tree}_trees.newick

${iqtree} -s ${ALI} -st DNA -m GTR+G -n 0 -z ${TREES} -zb 10000 -au -zw -nt 12

