#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 16
#$ -binding linear:16
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

DATASET="tara"
BACKBONE="Euks"
TAXCONSTR="Unconstr"
PACKAGE="${DATASET}_${BACKBONE}_${TAXCONSTR}_backbone"

source /etc/profile.d/modules.sh

module unload gcc
module load gcc/4.9.2
module load boost/1.60
module load openmpi/gcc

BASEDIR=/path/to/here

NUM_TASKS=16
export OMP_NUM_THREADS=16
export MP_SINGLE_THREAD=no
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

SAMPLEDIR="${BASEDIR}/06_samples/${PACKAGE}/samples_bplace"
WORKDIR="${BASEDIR}/08_emd/${PACKAGE}"

echo "Sample dir is ${SAMPLEDIR}"
echo "Work dir is   ${WORKDIR}"
echo "Running on `hostname`"
echo "`date` Start"
echo

cd ${BASEDIR}/software/genesis/bin
./bplace_emd ${SAMPLEDIR} ${WORKDIR}

echo
echo "`date` End"
