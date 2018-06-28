#!/bin/bash

genesis/bin/apps/silva_single_seqs Archaea &
PID_Archaea=$!

genesis/bin/apps/silva_single_seqs Bacteria &
PID_Bacteria=$!

genesis/bin/apps/silva_single_seqs Eukaryota &
PID_Eukaryota=$!

genesis/bin/apps/silva_single_seqs General &
PID_General=$!

wait ${PID_Archaea}
wait ${PID_Bacteria}
wait ${PID_Eukaryota}
wait ${PID_General}

BASEDIR="11_single_seqs/00_reference"

# genesis/bin/apps/fasta_cleanup ${BASEDIR}/Archaea_sequences.fasta
# genesis/bin/apps/fasta_cleanup ${BASEDIR}/Bacteria_sequences.fasta
# genesis/bin/apps/fasta_cleanup ${BASEDIR}/Eukaryota_sequences.fasta
# genesis/bin/apps/fasta_cleanup ${BASEDIR}/General_sequences.fasta

cd ${BASEDIR}
raxml -f c -m GTRGAMMA -s Archaea_sequences.cleaned.fasta -n Archaea_reduce >> Archaea_reduce.log
raxml -f c -m GTRGAMMA -s Bacteria_sequences.cleaned.fasta -n Bacteria_reduce >> Bacteria_reduce.log
raxml -f c -m GTRGAMMA -s Eukaryota_sequences.cleaned.fasta -n Eukaryota_reduce >> Eukaryota_reduce.log
raxml -f c -m GTRGAMMA -s General_sequences.cleaned.fasta -n General_reduce >> General_reduce.log
