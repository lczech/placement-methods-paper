#!/bin/bash

# REPLACE_METHOD 
# REPLACE_REFERENCE 
# REPLACE_CONSTR

rm -f top_eval.txt
rm -f top_eval_Constr.txt
rm -f top_eval_Unconstr.txt

DATADIR="../04_epa"

for METHOD in cavener majorities threshold_0.5 threshold_0.6 threshold_0.7 threshold_0.8 threshold_0.9 threshold_0.95 ; do
		
	echo "${METHOD}" >> top_eval.txt
	echo "===================" >> top_eval.txt
	echo >> top_eval.txt
			
	echo "${METHOD}" >> top_eval_Constr.txt
	echo "===================" >> top_eval_Constr.txt
	echo >> top_eval_Constr.txt	
		
	echo "${METHOD}" >> top_eval_Unconstr.txt
	echo "===================" >> top_eval_Unconstr.txt
	echo >> top_eval_Unconstr.txt
		
	for REFERENCE in Actinobacteria Cyanobacteria Proteobacteria Firmicutes Bacteroidetes ; do
		for CONSTR in Constr Unconstr ; do
		
			if [ -f "${DATADIR}/${METHOD}/${REFERENCE}_${CONSTR}/silva_tree_eval_tax_no-blacklist.log" ] ; then
				grep "Correct pqueries 0" ${DATADIR}/${METHOD}/${REFERENCE}_${CONSTR}/silva_tree_eval_tax_no-blacklist.log >> top_eval.txt
			else			
				grep "Correct pqueries 0" ${DATADIR}/${METHOD}/${REFERENCE}_${CONSTR}/silva_tree_eval_no-tax_no-blacklist.log >> top_eval.txt
			fi
			
			# same again, just split by constraint
			if [ -f "${DATADIR}/${METHOD}/${REFERENCE}_${CONSTR}/silva_tree_eval_tax_no-blacklist.log" ] ; then
				grep "Correct pqueries 0" ${DATADIR}/${METHOD}/${REFERENCE}_${CONSTR}/silva_tree_eval_tax_no-blacklist.log >> top_eval_${CONSTR}.txt
			else			
				grep "Correct pqueries 0" ${DATADIR}/${METHOD}/${REFERENCE}_${CONSTR}/silva_tree_eval_no-tax_no-blacklist.log >> top_eval_${CONSTR}.txt
			fi
			
		done
	done
	
	echo >> top_eval.txt
	echo >> top_eval_Constr.txt
	echo >> top_eval_Unconstr.txt
done
	
