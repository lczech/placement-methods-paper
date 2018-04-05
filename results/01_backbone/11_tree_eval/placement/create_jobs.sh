#!/bin/bash

rm -f submit.sh
echo "#!/bin/bash" > submit.sh

for CHUNK in `seq 0 5`; do
    mkdir job_${CHUNK}
    cat job_template.sh                            \
        | sed s/REPLACE_CHUNK/${CHUNK}/g           \
        > job_${CHUNK}/job_script.sh
    chmod 755 job_${CHUNK}/job_script.sh

    echo "cd job_${CHUNK}"    >> submit.sh
    echo "qsub job_script.sh" >> submit.sh
    echo "cd .."              >> submit.sh
    echo ""                   >> submit.sh
done

chmod 755 submit.sh
