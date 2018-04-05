#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich12 12
#$ -binding linear:12
#$ -q sandy.q
#$ -cwd
#$ -j y
#$ -l h_rt=08:00:00

echo "Running on `hostname`"
echo "Start at `date`"
echo

find * -maxdepth 0 -type d | xargs -I fn tar -czvf fn.tar.gz fn

echo
echo "End at `date`"
