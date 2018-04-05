#!/bin/bash

# Delete previous runs:
# rm -v */histogram*.csv
# rm -v */silva_tree_eval*.log

echo `date`
echo

#for METHOD in cavener majorities threshold_0.5 threshold_0.6 threshold_0.7 threshold_0.8 threshold_0.9 threshold_0.95 ; do
for METHOD in majorities ; do
	echo "==================================="
  	echo "    ${METHOD}"
 	echo "==================================="
 	
 	BASEDIR="11_consensus_seqs/04_epa/${METHOD}"
 	OUTDIR="11_consensus_seqs/05_accuracy_viz/${METHOD}"
 	mkdir -p ${OUTDIR}
 	
 	for REFERENCE in Archaea Bacteria Eukaryota General ; do
 		for CONSTR in Constr Unconstr ; do
 		
 			mkdir -p ${OUTDIR}/${REFERENCE}_${CONSTR}/
 			
 			#echo "Calling fuckcomma ${BASEDIR}/${REFERENCE}_${CONSTR}/epa_result.jplace"
 			#genesis/bin/apps/fuckcomma ${BASEDIR}/${REFERENCE}_${CONSTR}/epa_result.jplace
 		
 			#echo "Calling placement_histograms ${BASEDIR}/${REFERENCE}_${CONSTR}/epa_result.jplace ${OUTDIR}/${REFERENCE}_${CONSTR}"
  			#genesis/bin/apps/placement_histograms "${BASEDIR}/${REFERENCE}_${CONSTR}/epa_result.jplace" "${OUTDIR}/${REFERENCE}_${CONSTR}"
  			#echo 
 		done
 	done
 	

  	echo "Calling silva_tree_eval ${BASEDIR} ${OUTDIR}"
  	genesis/bin/apps/silva_tree_eval ${BASEDIR} ${OUTDIR}
  	echo
done


