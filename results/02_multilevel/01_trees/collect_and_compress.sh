#!/bin/bash

for METHOD in cavener majorities threshold_0.5 threshold_0.6 threshold_0.7 threshold_0.8 threshold_0.9 threshold_0.95 ; do

	for REFERENCE in Actinobacteria Cyanobacteria Proteobacteria Firmicutes Bacteroidetes ; do
	
		if [ ! -d "${METHOD}/${REFERENCE}_Unconstr" ]
		then
		  continue
		fi
		
		if [ -f "${METHOD}/${REFERENCE}_Unconstr/best_tree.newick" ]
		then
		  continue
		fi
	
		##############################
		#    Unconstrained
		##############################
		
		echo "======================================="
		echo "Processing ${METHOD}/${REFERENCE}_Unconstr"
		
		cd ${METHOD}/${REFERENCE}_Unconstr
		./collect_trees.sh
		tar -czvf seeds.tar.gz s[0-9]*
		cd ../..
		
	done
done

