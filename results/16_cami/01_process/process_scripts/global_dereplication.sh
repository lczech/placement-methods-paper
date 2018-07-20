#!/bin/bash

BASEDIR=/data/czechls/tara_amp
cd ${BASEDIR}

# Define binaries
VSEARCH=${BASEDIR}/software/vsearch-2.6.2-linux-x86_64/bin/vsearch
CUTADAPT="cutadapt"

# Define temporary files and output files
TMP_FASTA=tmp_fasta
FINAL_FASTA="all_derep.fasta"

# Pool sequences
cat fasta/*.fasta > "${TMP_FASTA}"

# Dereplicate (vsearch)
"${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
             --sizein \
             --sizeout \
             --fasta_width 0 \
             --output "${FINAL_FASTA}" > /dev/null

rm -f "${TMP_FASTA}"

