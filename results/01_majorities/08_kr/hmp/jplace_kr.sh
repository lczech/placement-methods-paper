#!/bin/bash

BASEDIR="/path/to/dir"

EPADIR=${BASEDIR}/06_samples/hmp/samples_jplace
OUTDIR=${BASEDIR}/08_kr/

echo "`date` Start"

${BASEDIR}/genesis/bin/apps/jplace_emd ${EPADIR} ${OUTDIR}

echo "`date` Finish"
