#!/bin/bash

set -e

name=${1:-cons}

echo $name

exit

OUT=profile/${name}/

mkdir -p ${OUT}
rm -fr ${OUT}*

for i in {0..63}; do
  echo "`date` Doing sample ${i}..."
  gappa analyze assign --jplace-path jplace/${name}/samples/sample_${i}.jplace \
    --taxon-file tree/bact_ref_pkg/tree_names_to_taxopath.txt \
    --taxonomy tree/bact_ref_pkg/ncbi_genera.tax --cami --sample-id ${i} \
    --out-dir ${OUT}/sample_${i} --verbosity 0


  cat ${OUT}/sample_${i}/cami.profile >> ${OUT}/merged.profile

done

echo "`date` All done!"
