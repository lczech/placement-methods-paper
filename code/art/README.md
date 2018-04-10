Overview
-------------------------

The programs in this directory were used for constructing 
the consensus sequences used for our Automatic Reference Trees method,
as well as for testing them. We used Silva version 123.1 for this.
The programs use genesis v0.19.0

`silva_consensus_seqs`
-------------------------

Our prototype implementation of `silva_entropy` (see below) only calculated
consensus sequences with the threshold method, using a 95% theshold.
In order to also use the other consensus methods, this program turns existing
consensus sequences into new ones, by reading the original alignment and
re-calculating sequences based on their names.

It would also have been possible to re-run the original program and just
replace the consensus method, but we wanted to make sure that we use the exact
same settings and sequences.

`silva_entropy`
-------------------------

Prototype implementation of the ART method.
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

`silva_tree_eval`
-------------------------

Main evaluation program of the ART trees, which measures the discrete and continuous
distances as reported in the paper, but also many other things.
Its output is a bit messy, but the most important files are written to the
directories `lists` and `tables`, which are then used by our plotting scripts.

`silva_tree_eval_blacklist`
-------------------------

Write a blacklist file that lists all sequence names (SEQ_xxxxxx_) 
that we want to try to exclude from evaluation.
This evaluation is mentioned in the ART supplement, where we wanted to assess
the effect of excluding "bad" sequences. We consider as bad:

 * all sequences form the Sativa mislabel list
 * all sequences that contain "incertae", "unclassified" or "unknown" in their taxopath

This list has 25,910 entries out of 598,470 sequences, or 4.3% of the data.
The results can be found in `data/silva/blacklist.txt`.

`silva_tree_eval_queries`
-------------------------

Preparation for the evaluation in `silva_tree_eval`.
The program writes out the ~600k silva sequences in to chunks, so that we can easily place 
them on the ARTs, and rename the sequences so that their names 
reflect the taxo path they should be placed on.

The renaming is done so that the evaulation is easy: simply look at the name to figure out where
it should go in the placement. This is because the names correspond to the taxonomic paths
of the tree leaves, and thus indicate where they should be placed (i.e., for which consensus
sequences they were used).
We also prepend a unique counter prefix to each sequence, so that
epa and other tools don't complain about duplicate sequence names (e.g., from the same genus).
This fixed length id is easy to remove in the eval program.
