Overview
-------------------------

This is the new run of the Bacteria backbone tree test with majority consensus sequences,
placing BV and HMP to get proper results for the plots in the paper.
The alignment and tree are copies from `11_consensus_seqs`,
where we tested the effects on accuarcy of using different consensus methods.
See there for how they were created.

`00_reference`: The reference alignment of the Bacteria from our ART method,
using majority rule consensus sequences.
It is only the Bact, as the other ones are not needed for this test.

`01_tree`: Reference tree inferred from the alignment.

`02_sequences`: The directory is not incldued here, but it just contained
links to the BV and HMP datasets, where the `chunks`, `filtered` and `maps` files are stored.

`03_align`: Aligning BV and HMP data to the reference tree and alignment.

`04_epa`: Placing the data to get jplace files.

`06_samples`: Similar to what we did in `01_backbone`, some files needed fixing
due to an issue with EPA-ng  the downside of being beta tester...
The directory contains those fixes, as well as the final unchunkified jplace files
for placing BV and HMP on the majority rule Bact tree.

`07_squash_edgepca`: Evaluation of the BV data using the two standard methods
Squash Clustering and Edge PCA. Those were used to evaluate how well the Bact tree
works for running standard placement analyses.

`08_kr`: KR distance calculations and evaluation on the HMP dataset,
again used to test how well the Bact tree works for real live data.
This contains the Edge PCA and MDS plots of the HMP data shown in the paper.
