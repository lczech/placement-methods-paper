#!/bin/bash
#
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -j y                              # join stout and stderr
#$ -q serial.q                       # Select the queue
#$ -l h_rt=24:00:00                  # Request the time for the job

. /etc/profile
. /etc/profile.d/modules.sh

BATCH=REPLACE_BATCH

total_lines=`cat ${BATCH} | wc -l`
line_counter=1

echo "Process batch ${BATCH}"
echo

for link in `cat ${BATCH}`; do
    file=${link##*/}

    echo "========================================="
    echo "    ${line_counter}/${total_lines}: ${file}"
    echo "========================================="

    cd ../fastq
    echo "Time: `date`"
    echo "Download: ${link}"
    wget -nv "ftp://${link}"
    #echo -n "Unpacking... "
    #gunzip ${file}
    #echo -e "done.\n"
    cd ../download_scripts

    line_counter=$(( $line_counter + 1 ))
done
