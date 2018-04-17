#!/bin/bash

BASEDIR="06_samples/tara_Gen_Unconstr_backbone"
GENESIS_BIN_DIR="software/genesis/bin/apps"

echo "Start at `date`"
echo

${GENESIS_BIN_DIR}/jplace_squash 6 ${BASEDIR}/samples_jplace ${BASEDIR}/squash

echo
echo "End at `date`"
