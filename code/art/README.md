Overview
=========================

The programs in this directory were used for constructing 
the consensus sequences used for our Automatic Reference Trees method,
as well as for testing them. We used Silva version 123 for this.
The programs use genesis v0.19.0

`silva_entropy.cpp`
-------------------------

This program reads the Silva sequences, calculates the entropy per taxon,
and builds consensus sequences. The parameters are hardcoded to the values
we used in the paper, but can be changed before compiling the program.

The resulting executable expects three input parameters:

 1. Taxonomy. In Silva 123, this is the file `tax_slv_ssu_123.1.txt`
 2. Sequence alignment. In Silva, this is the file 
    `SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta`,
    which has a width of 50.000 sites.
 3. An output dir to which to write the resulting files to.
 
The most important output of the program is the file `tax_cons_border.fasta`,
which contains the consensus sequences build by the algorithm.
The program also writes several other useful files with additional information.

`silva_consensus_seqs.cpp`
-------------------------

Our prototype implementation of `silva_entropy.cpp` (see above) only calculated
consensus sequences with the threshold method, using a 95% theshold.
In order to also use the other consensus methods, this program turns existing
consensus sequences into new ones, by reading the original alignment and
re-calculating sequences based on their names.

It would also have been possible to re-run the original program and just
replace the consensus method, but we wanted to make sure that we use the exact
same settings and sequences.
