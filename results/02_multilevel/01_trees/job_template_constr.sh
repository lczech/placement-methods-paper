#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 16
#$ -binding linear:16
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

module load gcc/4.9
module load boost/1.56_gcc


BASEDIR=/path/to/here

SATIVA_HOME=${BASEDIR}/software/sativa-master

ALI=${BASEDIR}/02_multilevel/00_reference/REPLACE_METHOD/Bacteria_REPLACE_REFERENCE.tips.reduced.fasta
TAX=${BASEDIR}/02_multilevel/00_reference/subtree_alignments/Bacteria_REPLACE_REFERENCE.tips.tax

${SATIVA_HOME}/epa_trainer.py -s ${ALI} -t ${TAX} -n REPLACE_REFERENCE_Constr_backbone -x ZOO -no-hmmer -N 20

#~ grep "\"raxmltree\":" Euks_Constr.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" > best_tree.newick
grep "\"raxmltree\":" REPLACE_REFERENCE_Constr_backbone.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" | sed "s/;\",[ ]*$/;/g" > best_tree.newick

