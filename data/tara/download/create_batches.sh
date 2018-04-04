#!/bin/bash

# basis: https://www.ebi.ac.uk/ena/data/view/PRJEB6610

mkdir -p ../fastq

NUM_JOBS=8

split -n l/${NUM_JOBS} "fastq_ftp_list" "batch_"

echo -e "#!/bin/bash\n" > submit.sh
chmod 755 submit.sh

for list in `ls batch_*`; do
    cat job_template.sh | sed s/REPLACE_BATCH/${list}/g > job_script_${list}.sh
    echo "qsub job_script_${list}.sh" >> submit.sh
done
