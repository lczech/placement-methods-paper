#!/bin/bash

BASEDIR=/data/czechls/tara_amp
VSEARCH=${BASEDIR}/software/vsearch-2.6.2-linux-x86_64/bin/vsearch
PEAR=${BASEDIR}/software/pear-0.9.11-linux-x86_64/bin/pear

THREADS="20"
MEM="1000M"

cd ${BASEDIR}
mkdir -p fastq_merged

for accession in `cat ${BASEDIR}/process_scripts/run_accessions.txt`; do
    echo "==============================="
    echo "    accession ${accession}"
    echo "==============================="
    echo "Got here at `date`"

    R1="${BASEDIR}/fastq_gz/${accession}_1.fastq.gz"
    R2="${BASEDIR}/fastq_gz/${accession}_2.fastq.gz"
    OUT_PREFIX="${BASEDIR}/fastq_merged/${accession}"
    
    # Wait for files to be downloaded
    LOADED=`grep "${accession}" ${BASEDIR}/download_scripts/download_finished.txt`
    while [ ! -n "${LOADED}" ]
	do
		sleep 2
		LOADED=`grep "${accession}" ${BASEDIR}/download_scripts/download_finished.txt`
	done
    #while [ ! -f ${R1} ]
	#do
	#	sleep 2
	#done
	#while [ ! -f ${R2} ]
	#do
	#	sleep 2
	#done
    
	# Last checks
    if [[ ! -e ${R1} ]] || [[ ! -e ${R2} ]] ; then
        echo "No fastq for ${accession}"
	    echo
        continue
    fi
    if [[ -e "${OUT_PREFIX}.assembled.fastq" ]] ; then
        echo "Merged file already exists for ${accession}"
	    echo
        continue
    fi
    
    echo "Started at `date`"

	echo "Decompress from gzip"
	zcat "${R1}" > "${R1/.gz/.raw}"

    # Test encoding
    echo "Test encoding"
    ENCODING=$(${VSEARCH} --fastq_chars "${R1/.gz/.raw}" 2>&1 | grep "^Guess" | grep -o "[0-9][0-9]$")
    if [[ "${ENCODING}" != 33 ]] && [[ "${ENCODING}" != 64 ]] ; then
        echo "Error: ${accession} unknown quality encoding: '${ENCODING}'" 1>&2
        continue
    fi
    
    echo "Delete decompressed file again"
    rm "${R1/.gz/.raw}"

	# Run pear
	echo "PEAR"
    ${PEAR} --forward-fastq ${R1} --reverse-fastq ${R2} --output ${OUT_PREFIX} --phred-base ${ENCODING} --memory ${MEM} --threads ${THREADS}
    #rm -f ${R1}
    #rm -f ${R2}
    
    echo "Compress merged files"
    cd ${BASEDIR}/fastq_merged
    tar -czf ${accession}.tar.gz ${accession}.*.fastq
    rm ${accession}.*.fastq
    cd ${BASEDIR}
    
    echo ${accession} >> process_scripts/pear_finished.txt

    echo "Finished at `date`"
    echo
done

