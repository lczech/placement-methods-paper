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
echo "`date` Starting EMD"
echo "================================================================="

cd ${BASEDIR}/software/genesis/bin/apps
./bplace_emd ${SAMPLEDIR} ${WORKDIR}
./mat_to_bmp ${WORKDIR}/emd.mat

###############################################################

echo "================================================================="
echo "`date` Starting EPCA"
echo "================================================================="

./bplace_epca ${SAMPLEDIR} ${WORKDIR}

###############################################################

echo "`date` Starting Viz"

echo "Each sample on its own"
for bplace in `ls ${SAMPLEDIR}/*.bplace`; do 
    name=${bplace#${SAMPLEDIR}/}
    name=${name%.bplace}
    echo
    echo "================================================================="
    echo "    ${name}"
    echo "================================================================="

    ./visualize_placements_svg ${bplace} ${WORKDIR}/samples/${name}
done

echo
echo "================================================================="
echo "    ALL"
echo "================================================================="

echo
echo "`date` All samples in one picture."
./visualize_placements_svg ${SAMPLEDIR} ${WORKDIR}/all

###############################################################

echo
echo "================================================================="
echo "    Zip"
echo "================================================================="

echo
echo "`date` Zipping sample pictures."

cd ${WORKDIR}
tar -czf samples.tar.gz samples/

echo
echo "`date` End"
