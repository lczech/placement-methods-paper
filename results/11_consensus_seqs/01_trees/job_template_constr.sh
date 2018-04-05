#!/bin/bash

module load gcc/4.9
module load boost/1.56_gcc


BASEDIR=/path/to/here


SATIVA_HOME=${BASEDIR}/software/sativa-master

ALI=${BASEDIR}/11_consensus_seqs/00_reference/REPLACE_METHOD/REPLACE_REFERENCE_sequences.fasta.reduced
TAX=${BASEDIR}/01_backbone/00_reference/REPLACE_REFERENCE/tax_assign.txt

${SATIVA_HOME}/epa_trainer.py -s ${ALI} -t ${TAX} -n REPLACE_REFERENCE_Constr_backbone -x ZOO -no-hmmer -N 20

#~ grep "\"raxmltree\":" Euks_Constr.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" > best_tree.newick
grep "\"raxmltree\":" REPLACE_REFERENCE_Constr_backbone.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" | sed "s/;\",[ ]*$/;/g" > best_tree.newick
