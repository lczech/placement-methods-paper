#!/bin/bash

module load gcc/4.9
module load boost/1.56_gcc

BASEDIR=/path/to/data
SATIVA_HOME=${BASEDIR}/software/sativa-master

ALI=${BASEDIR}/00_reference/General/Ref.phylip
TAX=${BASEDIR}/00_reference/General/tax_assign.txt

${SATIVA_HOME}/epa_trainer.py -s ${ALI} -t ${TAX} -n Gen_Constr_backbone -x ZOO -no-hmmer -N 40

#~ grep "\"raxmltree\":" Gen_Constr.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" > best_tree.newick
grep "\"raxmltree\":" Gen_Constr_backbone.refjson | sed "s/[ ]\+\"raxmltree\"\:\ \"//g" | sed "s/(r_/(/g" | sed "s/,r_/,/g" | sed "s/;\",[ ]*$/;/g" > best_tree.newick
