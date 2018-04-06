#!/bin/bash

for filename in `ls orig_queries/sample_*/*.jplace`; do
    # ls ${filename}
    # orig_queries/sample_p4z1r79/RAxML_portableTree.orig_queries.jplace

    prefix="orig_queries/sample_"
    suffix="/RAxML_portableTree.orig_queries.jplace"
    sample=${filename#$prefix}
    sample=${sample%$suffix}

    # echo ${sample}
    cp ${filename} orig_queries_jplace/${sample}.jplace
done
