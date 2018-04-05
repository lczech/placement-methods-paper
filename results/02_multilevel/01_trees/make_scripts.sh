#!/bin/bash

for METHOD in cavener majorities threshold_0.5 threshold_0.6 threshold_0.7 threshold_0.8 threshold_0.9 threshold_0.95 ; do

	mkdir -p ${METHOD}
	cp submit_jobs.sh ${METHOD}/submit_jobs.sh
	chmod 755 ${METHOD}/submit_jobs.sh

	for REFERENCE in Actinobacteria Cyanobacteria Proteobacteria Firmicutes Bacteroidetes ; do
	
		##############################
		#    Unconstrained
		##############################
		
		mkdir -p ${METHOD}/${REFERENCE}_Unconstr
		cat job_template_unconstr.sh | sed "s/REPLACE_METHOD/${METHOD}/g" | sed "s/REPLACE_REFERENCE/${REFERENCE}/g" > ${METHOD}/${REFERENCE}_Unconstr/job_template.sh
		
		cp collect_trees.sh ${METHOD}/${REFERENCE}_Unconstr/collect_trees.sh
		chmod 755 ${METHOD}/${REFERENCE}_Unconstr/collect_trees.sh
		
		##############################
		#    Constrained
		##############################

		mkdir -p ${METHOD}/${REFERENCE}_Constr
		
		cat job_template_constr.sh | sed "s/REPLACE_METHOD/${METHOD}/g" | sed "s/REPLACE_REFERENCE/${REFERENCE}/g" > ${METHOD}/${REFERENCE}_Constr/job_script.sh
		chmod 755 ${METHOD}/${REFERENCE}_Constr/job_script.sh

	done
done

