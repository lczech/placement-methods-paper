#!/bin/bash

# Purpose: scan both hmp data directories and produce uniq links to them in `links`.
# Leave out files that are too big (computation too long). 
# There are only 10 files >70MB, out of ~10,000 samples.
# It is not worth the effort to keep those two files in the analysis.

rm -rf links
mkdir -p links

for fn in `ls by_sample/`; do
    bn=${fn##*/}
    bn=${bn%.*}
    ex=${fn##*.}

    path="../by_sample/${fn}"

    #if [ `stat ${path} --printf="%s"` -lt "1000000" ]; then
    #    continue
    #fi
    if [ `stat by_sample/${fn} --printf="%s"` -gt "70000000" ]; then
        echo "Skip by_sample/${fn} (`stat by_sample/${fn} --printf="%s"` bytes)"
        continue
    fi
    ln -s ${path} links/smp_${bn}.${ex}
done

for fn in `ls SRP002860/`; do
    bn=${fn##*/}
    bn=${bn%.*}
    ex=${fn##*.}

    path="../SRP002860/${fn}"

    #if [ `stat ${path} --printf="%s"` -lt "1000000" ]; then
    #    continue
    #fi
    if [ `stat SRP002860/${fn} --printf="%s"` -gt "70000000" ]; then
        echo "Skip SRP002860/${fn} (`stat SRP002860/${fn} --printf="%s"` bytes)"
        continue
    fi
    ln -s ${path} links/srp_${bn}.${ex}
done
