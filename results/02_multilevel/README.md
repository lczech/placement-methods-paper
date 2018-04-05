Overview
-------------------------

Testing the accuracy of the second level placement, using a selection of five bacterial clades.
This is using different consensus methods, which was not reported in the paper. 
There, we just reported the results of using majority consensus for this experiment,
as this is what we also used for the other evals.
We only keep the majority rule results here to save space.

`00_reference`: Alignments for the five subclades, 
created using the program `code/art/silva_subtree_alignments.cpp`.

`01_trees`: Trees inferred from the alignmenets.

`02_sequences`: The selection of sequences for the five subclades,
created with the program `code/art/silva_subtree_sequences.cpp`.
We do not store the sequences here, as this is too much data.

`04_epa`: Placing the sequences on the trees.

`05_accuracy_viz`: Graph plotting of the accuracy.

`06_mean_genera`: As mentioned in the supplement, we also tested the influence 
of "mean" genera, which are hard to place, according to http://jcm.asm.org/content/45/9/2761.short
