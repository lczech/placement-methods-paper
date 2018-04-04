#!/bin/bash
#
#$ -S /bin/bash
#$ -pe mvapich12 12
#$ -binding linear:12
#$ -q sandy.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

DATASET="tara"
BACKBONE="Euks"
TAXCONSTR="Unconstr"
PACKAGE="${DATASET}_${BACKBONE}_${TAXCONSTR}_backbone"

BASEDIR="path/to/data"
MAP_DIR="${BASEDIR}/02_sequences/${DATASET}/maps/"
EPA_DIR="${BASEDIR}/06_samples/${PACKAGE}/chunks/"
OUT_DIR="${BASEDIR}/06_samples/${PACKAGE}/"

source /etc/profile.d/modules.sh
module unload gcc
module load gcc/4.9.3

# Copy and extract jplace files
echo "`date` Copy jplace files"
cd ${BASEDIR}/06_samples/${PACKAGE}/
cp -r ${BASEDIR}/04_epa/${PACKAGE}/jobs/ chunks
cd chunks
for f in *.tar.gz; do 
  tar xf $f
done
rm -f *.tar.gz
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

# Unchunkify
echo "`date` Unchunkify"
cd ~/genesis/bin
/usr/bin/time -v ./jplace_unchunkify ${MAP_DIR} ${EPA_DIR} ${OUT_DIR}
echo

# Count remaining pqueries
echo "`date` Remaining pqueries count in chunk jplace files:"
cd ${EPA_DIR}
rem_pqry_cnt=`find . -name "*.jplace" | xargs egrep "\"[0-9a-f]*\"" | wc -l`
echo

# Some more checks
echo "`date` Some more checks"
cd ${OUT_DIR}
echo "overview.csv length:   `wc -l overview.csv`"
echo "chunknames.csv length: `wc -l chunknames.csv`"
echo "missing.csv length:    `wc -l missing.csv`"
missing_uniq_cnt=`cat missing.csv | cut -f 1 | sort | uniq > missing.uniq`
echo "missing.uniq length:   `wc -l missing.uniq`"

echo "Remaining pqueries count plus missing.uniq length should equal number of initial reads. Please check manually:"
num_reads=$(( ${rem_pqry_cnt} + ${missing_uniq_cnt} ))
echo "${rem_pqry_cnt} + ${missing_uniq_cnt} = ${num_reads}"
echo

echo "Overview multiplicities: `awk -F"\t" '{sum+=$3} END{print sum;}' overview.csv`"
echo "Missing  multiplicities: `awk -F"\t" '{sum+=$2} END{print sum;}' missing.csv`"
echo

echo "Broken pquery list (obtained via regex) and missing.uniq (obtained from samples) should be the same. Using diff to check:"
echo "vvv"
diff chunks/matches.uniq missing.uniq
echo "^^^ empty?"
echo

echo "`date` Checking sample contents..."
mkdir check
for jf in `ls samples/*.jplace`; do
    sample_name=${jf#samples/}
    sample_name=${sample_name%.jplace}
    out_name="check/${sample_name}"

    echo "Checking sample ${sample_name}"

    # Build list of pquery names in newly created jplace file
    cat ${jf} | grep -o "\"[0-9a-f][0-9a-f]*\"" | grep -o "[0-9a-f]*" | sort > ${out_name}.smp

    # Build list of names from map file
    cat ${MAP_DIR}/${sample_name}.*.map | cut -f 1 | sort > ${out_name}.map

    # Compare: find the sequence names that are in the map, but not in the sample,
    # i.e., missing sequence names.
    diff ${out_name}.smp ${out_name}.map | egrep "[<>]" | grep -o "[0-9a-f]*" > ${out_name}.diff
done
echo "`date` done checking"
echo

# Now check the output of the above: missing sequnence names should be the ones that we removed previously
cat check/*.diff | sort | uniq > check.uniq

echo "The sequences missing in the samples (as compared to the original map files) should be the missing sequences that we previously removed. Thus, this should be empty:"
echo "vvv"
diff check.uniq missing.uniq
echo "^^^ empty?"
echo

# Move samples
echo "`date` Move samples"
mkdir samples_jplace
mkdir samples_bplace
mv samples/*.jplace samples_jplace
mv samples/*.bplace samples_bplace
chmod 444 samples_bplace/*
ls samples/
rm -r samples/
echo

# Compress and delete
echo "`date` Compress and delete unneeded stuff..."
tar -czf samples_jplace.tar.gz samples_jplace/
rm -r samples_jplace/

# tar -czf check.tar.gz check/
rm -r check/

tar -czf missing.tar.gz missing/
rm -r missing/

cd chunks/
find * -maxdepth 0 -type d | xargs -I fn tar -czf fn.tar.gz fn
ls -d */ | xargs rm -r
chmod 444 *

echo "`date` Done"
