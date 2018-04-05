#!/bin/bash

for i in 50
do
    mkdir run_${i}

    cat job_template.sh | sed s/REPLACE_NUM_RUNS/${i}/g > run_${i}/job_script_run_${i}.sh
    chmod 755 run_${i}/job_script_run_${i}.sh

    cd run_${i}
    
    # Submission accoring to cluster env
    qsub job_script_run_${i}.sh
    llsubmit job_script_run_${i}.sh
    cd ..
done
