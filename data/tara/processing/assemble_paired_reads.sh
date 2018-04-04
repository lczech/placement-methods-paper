#!/bin/bash

VSEARCH=vsearch
PEAR=pear

THREADS="32"
MEM="50M"

BATCH="run_accessions"
if [ "$#" -eq 1 ]; then
    BATCH=$1
fi

echo "Running batch ${BATCH}"
echo

BASEDIR="data/tara"
cd ${BASEDIR}

mkdir -p fastq_merged

for accession in `cat ${BATCH}`; do
    echo "==============================="
    echo "    accession ${accession}"
    echo "==============================="
    echo "start at `date`"

    R1="${BASEDIR}/fastq_gz/${accession}_1.fastq.gz"
    R2="${BASEDIR}/fastq_gz/${accession}_2.fastq.gz"
    OUT_PREFIX="${BASEDIR}/fastq_merged/${accession}"

    if [[ ! -e ${R1} ]] || [[ ! -e ${R2} ]] ; then
        echo "No fastq for ${accession}"
        continue
    fi

    if [[ -e "${OUT_PREFIX}.assembled.fastq" ]] ; then
        echo "Merged file already exists for ${accession}"
        continue
    fi

    # Test encoding
    ENCODING=$(${VSEARCH} --fastq_chars "${R1}" 2>&1 | grep "^Guess" | grep -o "[0-9][0-9]$")
    if [[ "${ENCODING}" != 33 ]] && [[ "${ENCODING}" != 64 ]] ; then
        echo "Error: ${accession} unknown quality encoding" 1>&2
        continue
    fi

    ${PEAR} --forward-fastq ${R1} --reverse-fastq ${R2} --output ${OUT_PREFIX} --phred-base ${ENCODING} --memory ${MEM} --threads ${THREADS}
    rm -f ${R1}
    rm -f ${R2}

    echo
done
