#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich16 320
#$ -binding linear:12
#$ -q bridge.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

module load gcc/4.9.3
module load boost/1.60
module load openmpi/gcc

BASEDIR=/path/to/data/
NUM_TASKS=12

BATCH=REPLACE_BATCH
WORKDIR="${BASEDIR}/01_majorities/04_epa/hmp/jobs/${BATCH}"

echo "This is the job for batch ${BATCH}"
echo "Workdir is ${WORKDIR}"
echo "Running on `hostname`"
echo "Start at `date`"
echo

EPA=${BASEDIR}/software/epa/bin/epa

# unpack alignment
echo
echo "Extracting alignments..."
cd ${BASEDIR}/01_majorities/03_align/hmp/jobs/
if [ ! -d "${BATCH}" ]; then
    tar xvzf ${BATCH}.tar.gz
fi
echo "Finished extracting"
echo

# Convert to fasta
echo "Converting alignments..."
for chunk in `ls ${BATCH}/papara_alignment.*`; do
    ${BASEDIR}/software/genesis/bin/phylip_fasta_conv ${BASEDIR}/01_majorities/03_align/hmp/jobs/${chunk}
done
echo "Finished converting"
echo

# Chunk list
cd ${BASEDIR}/01_majorities/03_align/hmp/jobs/${BATCH}
ls papara_alignment.*.fasta > chunk_list
QUERYDIR=${BASEDIR}/01_majorities/03_align/hmp/jobs/${BATCH}
#QUERIES=`awk -vORS=, '{ print "${QUERYDIR}/"$1 }' chunk_list | sed 's/,$/\n/'`
QUERIES=`awk -v QUERYDIR="${QUERYDIR}" -vORS=, '{ print QUERYDIR"/"$1 }' chunk_list | sed 's/,$/\n/'`
echo "queries: ${QUERIES}"
echo

# Input files
TREE=${BASEDIR}/01_majorities/01_tree/best_tree.newick
# ALI=${BASEDIR}/03_align/${PACKAGE}/jobs/${BATCH}/papara_alignment.${CHUNK}.fasta

# Run Forrest, run!
echo "Running EPA"

hostfile=${WORKDIR}/hostfile
awk '{printf("%s\n",$1);}' $PE_HOSTFILE > ${hostfile}
mpiexec --hostfile ${hostfile} -n 20 ${EPA} -g 0.98 -t ${TREE} -s ${QUERIES} -O -w ${WORKDIR}

echo "Finished EPA"

echo
echo "Files in align dir ${BASEDIR}/01_majorities/03_align/hmp/jobs/${BATCH}"
ls ${BASEDIR}/01_majorities/03_align/hmp/jobs/${BATCH}

echo
echo "Files in workdir ${WORKDIR}"
ls ${WORKDIR}

# deleting alignment
echo
echo "Deleting alignment..."
cd ${BASEDIR}/01_majorities/03_align/hmp/jobs/
#rm -r ${BATCH}/
echo "Finished deleting"

echo
echo "Batch job finished."
echo "End at `date`"
