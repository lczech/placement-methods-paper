Overview
-------------------------

The directory contains scripts and some result files for reproducing our
analyses from the papers. We do not include all result files here,
for storage reasons - compressed, those are more than 500GB of data...
Furthermore, we leave out many intermediate files, as even with the big ones
being compressed, we have 60k files in the project directory.
As per archiving policy, we have all those files backed up,
just in case those are needed in the future.

Thus, what you see here is a selection of the most imporant programs, scripts
and results, which should enable a reproduction of our results as reported
in the papers.

 * `01_backbone`: Most parts of the backbone/first level ARTs
   of the four trees of the paper (General, Archaea, Bacteria, Eukaryota).
   These however used a threshold consensus method,
   which is mostly good for the anlayses, but not for the final accuracy tests.
 * `01_majorities`: The final analyses of the BV and HMP datasets 
   using majority rule consensus sequences for creating the plots of the paper.
   This is just using the Bacteria tree.
 * `02_bact_clade_constraint`: Testing how the accuracy of the ART changes 
   when we do a high level tax constraint to get our five test clades constrained.
 * `02_multilevel`: Testing the accuracy of the second level placement,
   using a selection of five bacterial clades.
 * `11_consensus_seqs`: Testing the influence of using different consensus sequence
   methods on the accuracy of the ART placement.
