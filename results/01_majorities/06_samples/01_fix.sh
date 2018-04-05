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
EPA_DIR="${BASEDIR}/04_epa/${DATASET}/chunk"
OUT_DIR="${BASEDIR}/06_samples/${DATASET}/"

#source /etc/profile.d/modules.sh
#module unload gcc
#module load gcc/4.9.3

# Copy and extract jplace files
echo "`date` Copy jplace files"
cd ${OUT_DIR}/
mkdir -p chunks
cp -r ${EPA_DIR}/ chunks
cd chunks
#for f in *.tar.gz; do 
  #tar xf $f
#done
#rm -f *.tar.gz
echo

# Find broken pqueries
echo "`date` Find broken pqueries"
find . -name "*.jplace" | xargs perl -0777 -ne 'while(m/    \{"p":\n      \n      \],\n    "n": \["[0-9a-f]*"\]\n    \},?\n/g){print "$&\n";}' >> matches
find . -name "*.jplace" | xargs perl -0777 -ne 'while(m/,\n    \{"p":\n      \n      \],\n    "n": \["[0-9a-f]*"\]\n    \}\n/g){print "$&\n";}' >> weirdos
cat matches | grep -o "[0-9a-f]*" | tee matches.seq | sort | tee matches.sort | uniq > matches.uniq
echo

# Fix broken pquries
echo "`date` Fix broken pquries"
find . -name "*.jplace" | xargs perl -i -0777 -pe 's/    \{"p":\n      \n      \],\n    "n": \["[0-9a-f]*"\]\n    \},\n//g'
find . -name "*.jplace" | xargs perl -i -0777 -pe 's/,\n    \{"p":\n      \n      \],\n    "n": \["[0-9a-f]*"\]\n    \}\n//g'
echo

echo "`date` Done"
