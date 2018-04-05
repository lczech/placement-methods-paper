#!/bin/bash

JOBS=20

for REFERENCE in Archaea Bacteria Eukaryota General ; do
	
	##############################
	#    Unconstrained
	##############################
	
	cd ${REFERENCE}_Unconstr
	for i in `seq ${JOBS}`
	do
		SEED=$RANDOM

		while [ -d "s$SEED" ]; do
		    echo "Directory s$SEED already exists, skipping..."
		    SEED=$RANDOM
		done

		mkdir s${SEED}

		cat job_template.sh | sed s/REPLACE_SEED/${SEED}/g > s${SEED}/job_script_s${SEED}.sh
		chmod 755 s${SEED}/job_script_s${SEED}.sh

		cd s${SEED}
        qsub job_script_s${SEED}.sh
		llsubmit job_script_s${SEED}.sh
		cd ..
	done
	cd ..
	
	##############################
	#    Constrained
	##############################
	
	cd ${REFERENCE}_Constr
	qsub job_script.sh
	llsubmit job_script.sh
	cd ..
	
done

