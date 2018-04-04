#!/bin/bash

# Helper script to remove unnecessary gaps in the alignments.

RAXML=raxml

for domain in Archaea Bacteria Eukaryota General; do
    echo "${domain}: tax_cons_border.fasta"
    ${RAXML} -f c -m GTRGAMMA -s ../${domain}/tax_cons_border.fasta -n ${domain}_reduce > /dev/null
    mv ../${domain}/tax_cons_border.fasta.reduced ../${domain}/tax_cons_border.phylip
    ln -s tax_cons_border.phylip ../${domain}/Ref.phylip
done

