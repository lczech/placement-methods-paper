#!/bin/bash

# abort if any command fails
set -e

VSEARCH=vsearch
PEAR=pear
CUTADAPT=cutadapt

SCRIPTS=../place/scripts

HMM=${SCRIPTS}/mitags_extraction_protocol-master/miTAGs_extraction_protocol/HMM3/

THREADS=40
CLEAN=1

for P in "$@";
do
    # get the location to work in from args, ensure it is a folder and that the data is there
    SAMPLE_DIR=$(realpath ${P})
    WD=${SAMPLE_DIR}/reads

    cd ${WD}

    LOG=logfile
    FQ=${WD}/anonymous_reads.fq
    FQ_GZ=${FQ}.gz
    if [ -d ${WD} -a -r ${FQ_GZ} ]
    then echo "`date` Starting on " ${WD}
    else echo "`date` ERROR: Not a dir or no data: " ${WD} && exit 1
    fi

    # gunzip the interleaved fastq (preserve original)
    echo "`date` --- Decompressing $(basename ${FQ_GZ})..."
    gunzip < ${FQ_GZ} > ${FQ}

    # split the interleaved file
    R1=forward.fq
    R2=reverse.fq
    echo "`date` --- Deinterleaving " $(basename ${FQ})
    ${SCRIPTS}/deinterleave_fastq.sh < ${FQ} ${R1} ${R2}
    if [ -n ${CLEAN} ]; then rm ${FQ}; fi

    # check encoding via vsearch
    echo "`date` --- Test encoding..."
    ENCODING=$(${VSEARCH} --fastq_chars "${R1/.gz/.raw}" 2>&1 | grep "^Guess" | grep -o "[0-9][0-9]$")
    if [[ "${ENCODING}" != 33 ]] && [[ "${ENCODING}" != 64 ]] ; then
        echo "ERROR: unknown quality encoding: '${ENCODING}'" 1>&2
        exit 1
    fi

    # get some identifying name for the fasta file
    DIRNAME=$(basename ${SAMPLE_DIR})
    arr=(${DIRNAME//_/ })
    NAME=${arr[2]}_${arr[3]}

    # Run pear
    echo "`date` --- Paired-end read merging..."
    OUT_PREFIX=merged
    ${PEAR} --forward-fastq ${R1} --reverse-fastq ${R2} --output ${OUT_PREFIX} --phred-base ${ENCODING} --threads ${THREADS} 2>> "${LOG}"
    if [ -n ${CLEAN} ]; then rm ${R1} ${R2}; fi

    ASSEMB=${OUT_PREFIX}.assembled.fastq
    FA=${NAME}.fasta

    echo "`date` --- Discard erroneous sequences and convert to fasta..."
    "${VSEARCH}" \
        --quiet \
        --fastq_filter "${ASSEMB}" \
        --fastq_maxns 0 \
        --fastaout "${FA}" 2>> "${LOG}"

    if [ -n ${CLEAN} ]; then rm ${OUT_PREFIX}.*; fi

    # index the fasta file
    echo "`date` --- Indexing..."
    cdbfasta ${FA}

    FA_RNA=${NAME}.rRNA
    # select ssu seqs
    echo "`date` --- HMMing..."
    ${SCRIPTS}/rna_hmm3.py -i ${FA} -o ${FA_RNA} -m ssu -k bac,arc,euk -L ${HMM} -p ${THREADS}

    echo "`date` --- Parsing..."
    ${SCRIPTS}/parse_rna_hmm3_output.pl ${FA_RNA}

    echo "`date` --- Extracting..."
    ${SCRIPTS}/extract_rrna_seqs.pl ${FA_RNA}.parsed 1 100

    echo "`date` Done!"

    cd -

done
