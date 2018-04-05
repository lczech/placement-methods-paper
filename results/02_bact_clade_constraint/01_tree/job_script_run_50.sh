#!/bin/bash
#
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -q sandy.q                        # Select the queue
#$ -pe mvapich12 12                  # Set the parallel environment
#$ -l h_rt=48:00:00                  # Request the time for the job


module load gcc/4.9
module load boost/1.56_gcc

BASEDIR=/path/to/here/

NUM_RUNS=50

SATIVA_HOME=${BASEDIR}/software/sativa-master

ALI=${BASEDIR}/02_bact_clade_constraint/00_reference/Bacteria_sequences.fasta.reduced
TAX=${BASEDIR}/02_bact_clade_constraint/00_reference/tax_assign_trimmed.txt

${SATIVA_HOME}/epa_trainer.py -s ${ALI} -t ${TAX} -n Bact_Clade_Constr_${NUM_RUNS} -x ZOO -no-hmmer -N ${NUM_RUNS}

#~ grep "\"raxmltree\":" Euks_Constr.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" > best_tree.newick
grep "\"raxmltree\":" Bact_Clade_Constr_${NUM_RUNS}.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" | sed "s/;\",[ ]*$/;/g" > best_tree_${NUM_RUNS}.newick
