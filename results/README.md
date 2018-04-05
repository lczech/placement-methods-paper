Overview
-------------------------

The directory contains scripts and some result files for reproducing our
analyses from the papers. We do not include all result files here,
for storage reasons - compressed, those are more than 500GB of data...

 * `01_backbone`: Most parts of the backbone/first level ARTs
   of the four trees of the paper (General, Archaea, Bacteria, Eukaryota).
   These however used a threshold consensus method,
   which is mostly good for the anlayses, but not for the final accuracy tests.
 * `01_majorities`: The final analyses of the BV and HMP datasets 
   using majority rule consensus sequences for creating the plots of the paper.
   This is just using the Bacteria tree.
 * `02_bact_clade_constraint`: Test how the accuracy of the ART changes 
   when we do a high level tax constraint to get our five test clades constrained.
