#!/bin/bash

# REPLACE_METHOD 
# REPLACE_REFERENCE 
# REPLACE_CONSTR

for METHOD in cavener majorities threshold_0.5 threshold_0.6 threshold_0.7 threshold_0.8 threshold_0.9 threshold_0.95 ; do
	for REFERENCE in Archaea Bacteria Eukaryota General ; do
		for CONSTR in Constr Unconstr ; do
		
			echo "======================================================"
			
			if [ ! -f "../01_trees/${METHOD}/${REFERENCE}_${CONSTR}/best_tree.newick" ] ; then
				echo "Skipping ${METHOD} ${REFERENCE} ${CONSTR}. No Tree: ../01_trees/${METHOD}/${REFERENCE}_${CONSTR}/best_tree.newick"
				continue
			fi
			if [ -f "${METHOD}/${REFERENCE}_${CONSTR}/job.sh" ] ; then
				echo "Skipping ${METHOD} ${REFERENCE} ${CONSTR}. Already done or running: ${METHOD}/${REFERENCE}_${CONSTR}/job.sh"
				continue
			fi
			
			mkdir -p ${METHOD}/${REFERENCE}_${CONSTR}
			cat job_template.sh | sed "s/REPLACE_METHOD/${METHOD}/g" | sed "s/REPLACE_REFERENCE/${REFERENCE}/g" | sed "s/REPLACE_CONSTR/${CONSTR}/g" > ${METHOD}/${REFERENCE}_${CONSTR}/job.sh
			chmod 755 ${METHOD}/${REFERENCE}_${CONSTR}/job.sh
			
			cd ${METHOD}/${REFERENCE}_${CONSTR}
			
			echo "Submitting ${METHOD}/${REFERENCE}_${CONSTR}/job.sh"
			llsubmit job.sh
			cd ../..
			
		done
	done
done
	
