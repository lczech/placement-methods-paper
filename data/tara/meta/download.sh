#!/bin/bash

# basis: https://www.ebi.ac.uk/ena/data/view/PRJEB6610

mkdir -p xml

total_lines=`cat list.txt | wc -l`
line_counter=1

while IFS='' read -r line || [[ -n "$line" ]]; do
    src=`echo $line | cut -d' ' -f 1`
    dst=`echo $line | cut -d' ' -f 3`

    echo -e "\e[32m========================================="
    echo "    ${line_counter}/${total_lines}: ${dst}"
    echo -e "=========================================\e[0m"

    cd xml
    echo "Time: `date`"
    echo "Download: ${src}"
    wget -O ${dst} "http://www.ebi.ac.uk/ena/data/view/${src}&display=xml&download=xml&filename=${src}.xml"
    #echo -n "Unpacking... "
    #gunzip ${file}
    #echo -e "done.\n"
    cd ..

    line_counter=$(( $line_counter + 1 ))
done < "list.txt"


