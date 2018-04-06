#!/bin/bash

DATASET="tara"
BACKBONE="Gen"
TAXCONSTR="Constr"
PACKAGE="${DATASET}_${BACKBONE}_${TAXCONSTR}_backbone"

BASEDIR=/path/to/here

NUM_TASKS=32
export OMP_NUM_THREADS=32
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

SAMPLEDIR="${BASEDIR}/06_samples/${PACKAGE}/samples_bplace"
WORKDIR="${BASEDIR}/08_emd/${PACKAGE}"

mkdir -p ${WORKDIR}/samples

echo "Sample dir is ${SAMPLEDIR}"
echo "Work dir is   ${WORKDIR}"
echo "Running on `hostname`"
echo "`date` Start"
echo

###############################################################

echo "================================================================="
echo "`date` Starting Matrices"
echo "================================================================="

cd ${BASEDIR}/software/genesis/bin/apps
./bplace_matrices ${SAMPLEDIR} ${WORKDIR}

###############################################################

echo "================================================================="
echo "`date` Starting NHD"
echo "================================================================="

./bplace_nhd ${SAMPLEDIR} ${WORKDIR}

###############################################################

echo
echo "`date` End"
