#!/bin/bash
#
#$ -S /bin/bash                      # Use bash
#$ -cwd -V                           # Shift directories and export variables
#$ -j y                              # join stout and stderr
#$ -q serial.q                       # Select the queue
#$ -l h_rt=24:00:00                  # Request the time for the job

. /etc/profile
. /etc/profile.d/modules.sh

cd data/tara/fastq_bz2/

for f in *.fastq.gz ; do
    printf "%s\n" "${f}"
    zcat "${f}" | bzip2 -z9fc > "${f/.gz/.bz2}"
    rm -f "${f}"
done
