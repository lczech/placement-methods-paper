#!/bin/bash

BASEDIR=/path/to/here

PLACE_FILE=${BASEDIR}/02_bact_clade_constraint/02_epa/epa_result.jplace
OUT_DIR=${BASEDIR}/02_bact_clade_constraint/03_eval
PROG=${BASEDIR}/genesis/bin/apps/silva_subclade_hits

${PROG} ${PLACE_FILE} ${OUT_DIR} &> eval.log
