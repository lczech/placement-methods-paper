#!/bin/bash

BASEDIR=`pwd`

for filename in `ls orig_queries_bplace/*.jplace`; do
    echo ${filename}

    cd genesis/bin/apps
    ./jplace_to_bplace ${BASEDIR}/${filename}
    cd ${BASEDIR}
done
