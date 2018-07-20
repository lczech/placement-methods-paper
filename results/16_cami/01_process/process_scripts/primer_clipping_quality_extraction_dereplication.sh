#!/bin/bash

BASEDIR=/data/czechls/tara_amp
cd ${BASEDIR}

# Create output dirs
mkdir -p fasta
#mkdir -p qual
mkdir -p fasta_log

# Define variables
PRIMER_F="TTGTACACACCGCCC"
PRIMER_R="GTAGGTGAACCTGCNGAAGG"
MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))  # primer match is >= 2/3 of primer length
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))
#QUALITY_FILE="qual/TARA_V9_samples.qual"

# Define binaries
VSEARCH=${BASEDIR}/software/vsearch-2.6.2-linux-x86_64/bin/vsearch
CUTADAPT="cutadapt"

# Define temporary files and output files
INPUT_REVCOMP=tmp_input_revcomp
TMP_FASTQ=tmp_fastq1
TMP_FASTQ2=tmp_fastq2
TMP_FASTA=tmp_fasta
OUTPUT=tmp_output

for accession in `cat ${BASEDIR}/process_scripts/run_accessions.txt`; do
    echo "==============================="
    echo "    accession ${accession}"
    echo "==============================="
    echo "Got here at `date`"

	# Wait for files to be peared
    LOADED=`grep "${accession}" ${BASEDIR}/process_scripts/pear_finished.txt`
    while [ ! -n "${LOADED}" ]
	do
		sleep 2
		LOADED=`grep "${accession}" ${BASEDIR}/process_scripts/pear_finished.txt`
	done

	INPUT_TAR="fastq_merged/${accession}.tar.gz"
    
    # Safety
    if [[ ! -e ${INPUT_TAR} ]] ; then
        echo "No fastq tar for ${accession}"
        echo
        continue
    fi
    if [[ -e "fasta/${accession}.fasta" ]] ; then
        echo "Fasta file already exists for ${accession}"
        echo
        continue
    fi
    
    echo "Started at `date`"
    
    echo "Uncompress peared fastq files"
    cd fastq_merged
    tar xf ${accession}.tar.gz
	cd ${BASEDIR}
	INPUT="fastq_merged/${accession}.assembled.fastq"
	
	# Last check
	if [[ ! -e ${INPUT} ]] ; then
        echo "No fastq for ${accession}"
        echo
        continue
    fi

    # Reverse complement fastq file
    echo "`date` Reverse complement fastq file"
    "${VSEARCH}" \
        --quiet \
        --fastx_revcomp "${INPUT}" \
        --fastqout "${INPUT_REVCOMP}"

    LOG="fasta_log/${accession}.log"
    FINAL_FASTA="fasta/${accession}.fasta"

    # Trim tags, forward & reverse primers (search normal and antisens)
    echo "`date` Trim tags, forward & reverse primers (search normal and antisens)"
    cat "${INPUT}" "${INPUT_REVCOMP}" | \
        ${CUTADAPT} --discard-untrimmed --minimum-length ${MIN_LENGTH} -g "${PRIMER_F}" -O "${MIN_F}" - 2>> "${LOG}" | \
        ${CUTADAPT} --discard-untrimmed --minimum-length ${MIN_LENGTH} -a "${PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${TMP_FASTQ}"

    # Discard erroneous sequences and add expected error rates
    echo "`date` Discard erroneous sequences and add expected error rates"
    "${VSEARCH}" \
        --quiet \
        --fastq_filter "${TMP_FASTQ}" \
        --fastq_maxns 0 \
        --relabel_sha1 \
        --eeout \
        --fastqout "${TMP_FASTQ2}" 2>> "${LOG}"

    # Convert fastq to fasta (discard sequences containing Ns)
    echo "`date` Convert fastq to fasta (discard sequences containing Ns)"
    "${VSEARCH}" \
        --quiet \
        --fastq_filter "${TMP_FASTQ}" \
        --fastq_maxns 0 \
        --fastaout "${TMP_FASTA}" 2>> "${LOG}"

    # Dereplicate at the study level (vsearch)
    echo "`date` Dereplicate at the study level (vsearch)"
    "${VSEARCH}" \
        --quiet \
        --derep_fulllength "${TMP_FASTA}" \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --output "${FINAL_FASTA}" 2>> "${LOG}"

    # Discard quality lines, extract sha1, expected error rates and read length
    echo "`date` Discard quality lines, extract sha1, expected error rates and read length"
    sed 'n;n;N;d' "${TMP_FASTQ2}" | \
        awk 'BEGIN {FS = "[;=]"}
             {if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
        tr -d "@" >> "${OUTPUT}"

    # Erase fastq file
    #rm "${INPUT}"
    rm fastq_merged/${accession}.*.fastq
    
    echo ${accession} >> process_scripts/fasta_finished.txt

	echo "Finished at `date`"
    echo
done

# Produce the final quality file
#sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
#    uniq --check-chars=40 > "${QUALITY_FILE}"

# Clean
rm -f "${INPUT_REVCOMP}" "${TMP_FASTQ}" "${TMP_FASTA}" "${TMP_FASTQ2}" "${OUTPUT}"

