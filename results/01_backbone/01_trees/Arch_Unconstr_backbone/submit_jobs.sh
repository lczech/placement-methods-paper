#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Please specify the number of jobs as the first argument!"
  exit
else
  JOBS=$1
fi


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
    
    # Submit script according to cluster settings
    #qsub job_script_s${SEED}.sh
    #llsubmit job_script_s${SEED}.sh
    cd ..
done
