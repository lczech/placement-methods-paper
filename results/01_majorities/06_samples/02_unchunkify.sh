#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich12 12
#$ -binding linear:12
#$ -q sandy.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

DATASET="bv"
PACKAGE="bv"

BASEDIR="/path/to/data/01_majorities"
MAP_DIR="${BASEDIR}/02_sequences/${DATASET}/maps/"
EPA_DIR="${BASEDIR}/06_samples/${DATASET}/chunks/"
OUT_DIR="${BASEDIR}/06_samples/${DATASET}/"

# Unchunkify
echo "`date` Unchunkify"
cd ${BASEDIR}/genesis/bin
/usr/bin/time -v ./jplace_unchunkify ${MAP_DIR} ${OUT_DIR} ${OUT_DIR}
echo

echo "`date` Done"
