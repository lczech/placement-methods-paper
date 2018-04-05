#!/bin/bash

BASEDIR="/path/to/dir"

EPCA=${BASEDIR}/01_majorities/genesis_new/bin/apps/jplace_epca_vis

EPA_DIR=${BASEDIR}/01_majorities/06_samples/hmp/samples_jplace
OUT_DIR=${BASEDIR}/01_majorities/08_kr/hmp

${EPCA} ${EPA_DIR} ${OUT_DIR}
