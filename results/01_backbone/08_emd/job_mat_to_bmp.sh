#!/bin/bash

BASEDIR=/path/to/here
WORKDIR=`pwd`

for file in `ls */*.mat`; do
    echo ${WORKDIR}/${file}
    ${BASEDIR}/software/genesis/bin/apps/mat_to_bmp ${WORKDIR}/${file}
done
