Overview
-------------------------

The programs in this directory are the prototypes of the multilevel methods.

`fasta_chunkify`
-------------------------

Prototype of the `chunkify` command in gappa.
It takes an input and an output path, and processes all fasta files in the input
by splitting them into chunks of 50k sequences, and writing abundance maps for each of them
Furthermore, because it is a dataset specific prototype,
it does some processing and filtering of the HMP and Tara data
that would otherwise have been extra functions.
For example, sequences out of our specified length requirements are exlcuded from the output.

`fasta_chunkify_simple`
-------------------------

Split a fasta file into chunks, but without abundance maps.

`jplace_unchunkify`
-------------------------

The program takes two command line arguments:

 1. Abundance map file directory
 2. Jplace directory with the chunified jplace files.
 3. Output directory to write files to.

It then undoes the chunkfiy process, that is, it writes per-sample jplace files
with their abundances per sequence.

`silva_subtree_alignments`
-------------------------

Write the alignments of consensus sequences and the taxonomies for the five clades
that we want to evaluate for the Russian doll approach.
Used in `results/02_multilevel`.

`silva_subtree_sequences`
-------------------------

Write fasta files that contain the aligned sequences of the five clades that
we want to evaluate for the Russian doll approach.

