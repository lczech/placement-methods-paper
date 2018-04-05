#!/bin/bash

# Delete previous runs:
# rm -v */histogram*.csv
# rm -v */silva_tree_eval*.log

echo `date`
echo

GENESIS_BIN="/path/to/genesis/bin/apps"
  
for dir in * ; do
  if [ ! -e "${dir}/epa_result.jplace" ]; then
    #echo -e "Skipping\n"
    continue
  fi

  echo "==================================="
  echo "    ${dir}"
  echo "==================================="

  taxon=`echo ${dir} | sed 's/_.*//'`
  if [ "${taxon}" == "Euks" ]; then
    # Unfortunate special case...
    taxon="Euk"
  fi
  if [ "${taxon}" == "Gen" ]; then
    # Unfortunate special case...
    taxon=""
  fi
  echo -e "Taxon: ${taxon}\n"

  if [ -e "${dir}/silva_tree_eval_no-tax_no-blacklist.log" ]; then
    echo -e "Already done\n"
    continue
  fi

  echo "Calling placement_histograms ${dir}/epa_result.jplace ${dir}"
  ${GENESIS_BIN}/placement_histograms "${dir}/epa_result.jplace" ${dir}
  echo 

  echo "Calling silva_tree_eval ${dir}/epa_result.jplace ${taxon}"
  ${GENESIS_BIN}/silva_tree_eval "${dir}/epa_result.jplace" ${taxon}
  echo
done
