#!/bin/bash

for BEST_TREE in `find . -name "best_tree.newick"` ; do
	
	BASE=`dirname ${BEST_TREE}`
	SIZE=`stat --printf="%s" ${BEST_TREE}`
	
	if [ "${SIZE}" -eq "0" ] ; then
		echo "Resubmitting ${BASE}" 
		
		rm -rf ${BASE}/best_tree.newick
		rm -rf ${BASE}/job*.err
		rm -rf ${BASE}/job*.out
		
		cd ${BASE}
		llsubmit job_script.sh
		cd - &> /dev/null
	fi
	
done

