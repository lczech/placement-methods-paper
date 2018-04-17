#!/bin/bash

BASEDIR="06_samples/bv_Gen_Unconstr_backbone"
GENESIS_BIN_DIR="software/genesis/bin/apps"

echo "Start at `date`"
echo

mkdir ${BASEDIR}/squash

${GENESIS_BIN_DIR}/bplace_squash 4 ${BASEDIR}/samples_bplace ${BASEDIR}/squash

cd ${BASEDIR}/squash
mkdir svg_c
mkdir svg_r
mkdir nexus

mv *.nexus nexus
mv *_c.svg svg_c
mv *_r.svg svg_r

find * -maxdepth 0 -type d | xargs -I fn tar -czf fn.tar.gz fn

echo
echo "End at `date`"
