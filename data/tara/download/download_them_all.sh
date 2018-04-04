#!/bin/bash

# basis: https://www.ebi.ac.uk/ena/data/view/PRJEB6610

mkdir -p fastq

BATCH=$1
echo -e "Batch ${BATCH}\n"

total_lines=`cat ${BATCH} | wc -l`
line_counter=1

for link in `cat ${BATCH}`; do
    file=${link##*/}

    echo -e "\e[32m========================================="
    echo "    ${line_counter}/${total_lines}: ${file}"
    echo -e "=========================================\e[0m"

    cd fastq
    echo "Time: `date`"
    echo "Download: ${link}"
    wget "ftp://${link}"
    #echo -n "Unpacking... "
    #gunzip ${file}
    #echo -e "done.\n"
    cd ..

    line_counter=$(( $line_counter + 1 ))
done
