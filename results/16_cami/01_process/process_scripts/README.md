after or while downloading, we can start the two step assemble scripts:
first, assemble_paired_reads.sh to run PEAR, 
then primer_clipping_quality_extraction_dereplication.sh to get the final fasta files.
in the end, we put the amplicon couns per sequence into a simpler format by running
`sed -i -- 's/;size=/_/g; s/;//g' *` in the `fasta` dir.
