#!/bin/bash

BASEDIR="/path/to/here"

echo "Start at `date`"
echo

for threads in 32 16 8 4 2 1 ; do
    echo "===================================================="
    echo "Running on ${thread} threads"

    /usr/bin/time -v ${BASEDIR}/genesis/bin/apps/jplace_emd_speed_comp ${threads} ${BASEDIR}/orig_queries_jplace_clean ${BASEDIR}/
    echo
done

echo
echo "End at `date`"
