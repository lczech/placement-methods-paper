#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 16
#$ -binding linear:16
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=08:00:00

DATASET="tara"
BACKBONE="Euks"
TAXCONSTR="Unconstr"
PACKAGE="${DATASET}_${BACKBONE}_${TAXCONSTR}_backbone"

module unload mpi.ibm/1.4
module unload gcc
module load gcc/4.9
module load boost/1.56_gcc
module load mpi.intel/5.1_gcc

BASEDIR=/path/to/here

NUM_TASKS=28
export OMP_NUM_THREADS=28
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

SAMPLEDIR="${BASEDIR}/06_samples/${PACKAGE}/samples_bplace"
WORKDIR="${BASEDIR}/09_viz/${PACKAGE}"

mkdir -p ${WORKDIR}/samples

echo "Sample dir is ${SAMPLEDIR}"
echo "Work dir is   ${WORKDIR}"
echo "Running on `hostname`"
echo "`date` Start"
echo

cd ${BASEDIR}/software/genesis/bin

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
echo "All samples in one picture."
./visualize_placements_svg ${SAMPLEDIR} ${WORKDIR}/all

echo
echo "`date` End"
