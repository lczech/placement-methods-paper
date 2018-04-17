Overview
-------------------------

This directory contains the resulting alignments from our ART method for the four 
trees `Archaea`, `Bacteria`, `Eukaryota` and `General`.

The files were created with the program `code/art/silva_entropy.cpp`.
See there for details on the processing and the input files.
For storage space reasons, we here only keep the most important files per run:

 * `tax_cons_border.*`: The main alignment file of the consensus 
   sequences for the taxonomic ranks that were selected by the algorithm.
 * `tax_assign.txt`: List of all selected taxa for which there are consensus sequences.
 * `tax_pruned.txt`: Full list of the selected taxa as well as all higher level ranks.
 
Furthermore, `tax_all_entr.csv` contains the entropy calculation of all taxa.
This can be used to check which ranks have which entropy.

Lastly `entropy` contains some tests of how the entropy measurement behaves
with our files, and as a measure in general
