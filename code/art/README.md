Overview
-------------------------

The programs in this directory were used for constructing 
the consensus sequences used for our Automatic Reference Trees method,
as well as for testing them. We used Silva version 123.1 for this.
The programs use genesis v0.19.0

`fasta_chunkify`
-------------------------

Prototype of the `chunkify` command in gappa.
It takes an input and an output path, and processes all fasta files in the input
by splitting them into chunks of 50k sequences, and writing abundance maps for each of them
Furthermore, because it is a dataset specific prototype,
it does some processing and filtering of the HMP and Tara data
that would otherwise have been extra functions.
For example, sequences out of our specified length requirements are exlcuded from the output.

`silva_consensus_seqs`
-------------------------

Our prototype implementation of `silva_entropy.cpp` (see above) only calculated
consensus sequences with the threshold method, using a 95% theshold.
In order to also use the other consensus methods, this program turns existing
consensus sequences into new ones, by reading the original alignment and
re-calculating sequences based on their names.

It would also have been possible to re-run the original program and just
replace the consensus method, but we wanted to make sure that we use the exact
same settings and sequences.

`silva_entropy`
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
The results of this program for our four test cases
`Archaea`, `Bacteria`, `Eukaryota` and `General` are stored in 
`results/01_backbone/00_reference`. See there for details.

`silva_subclade_hits`
-------------------------

Used in `results/02_bact_clade_constraint/03_eval` to test how many sequences
were placed in their correct clade when using the clade-constrained tree.

`silva_subtree_alignments`
-------------------------

Write the alignments of consensus sequences and the taxonomies for the five clades
that we want to evaluate for the Russian doll approach.
Used in `results/02_multilevel`.

`silva_subtree_sequences`
-------------------------

Write fasta files that contain the aligned sequences of the five clades that
we want to evaluate for the Russian doll approach.
