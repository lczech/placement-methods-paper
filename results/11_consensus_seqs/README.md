Overview
-------------------------

Here, we test the accuarcy of the ART trees using different consenus methods:

 * cavener
 * majorities
 * threshold_0.5
 * threshold_0.6
 * threshold_0.7
 * threshold_0.8
 * threshold_0.9
 * threshold_0.95

Trees are inferred contrained and unconstrained.
That is, a lot of combinations are tested, so we cannot publish all results here.
In the paper, we only reported some of the methods,
and excluded many of the thresold values, as well as the constraind variants of all of them.
Here, you find the full list, and all figures.

`00_reference`: Reference sequences for the four trees General, Archaea, Bact and Euks,
for each of the tested consensus sequence methods. These alignments were created with `code/art/silva_consensus_seqs.cpp`.
We only keep the scripts here, because of storage...

`01_trees`: Infer trees for all combinations.

`04_epa`: Place the silva sequences on all the trees, to test their accuracy.

`05_accuracy_viz`: Visualization, that is, the plots for this test as found in the paper,
as well as all the ones that we left out (constrained trees and other threshold methods).

`06_subclades`: Testing if all the five subclades that we used for the multilevel tests are good.
That is, we see how many of the placed sequences ended up in their respective correct clade,
using the program `code/art/silva_subclade_hits 11_consensus_seqs/04_epa/majorities/Bacteria_Unconstr/epa_result.jplace 11_consensus_seqs/06_subclades/ > eval.log`

The resulting accuracy:

Bacteria_Actinobacteria_:	51035	/	51160	=	0.997557
Bacteria_Firmicutes_:	131213	/	138517	=	0.94727
Bacteria_Proteobacteria_:	199121	/	200083	=	0.995192
Bacteria_Bacteroidetes_:	48788	/	49174	=	0.99215
Bacteria_Cyanobacteria_:	11316	/	11379	=	0.994463
total	441473		/	450313		0.9803692099

Furthermore, the directory contains the clade tree figure of the paper.

`13_rf_distances`: Test to get the rf distances between trees of difference consensus methods.
