#!/bin/bash

BASEDIR=bacterial_vaginosis

for sample_name in `ls ${BASEDIR}/00_data/queries/`
do
  sample_name=$(echo ${sample_name} | cut -f 1 -d '.')
  mkdir sample_${sample_name}

  cat job_template.sh | sed s/REPLACE_SAMPLE_NAME/${sample_name}/g > sample_${sample_name}/job_script.sh

  echo "Submitting job ${sample_name}..."

  # submit depending on cluster env
  cd sample_${sample_name}
  qsub job_script.sh
  llsubmit job_script.sh
  cd ..

done
