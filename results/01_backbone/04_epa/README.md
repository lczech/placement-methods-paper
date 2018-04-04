Overview
-------------------------

This directory contains the scripts of running EPA on the query sequences 
with the reference consensus sequences and the reference tree inferred from them.
We use EPA-ng for this step.
Furthermore, we use a genesis program called `phylip_fasta_conv`,
which is now located in `code/tools/phylip_to_fasta.cpp` in order to convert 
the papara output Phylip files to fasta.

We do not offer all placement result files here, as this would be too much data.
Instead, we provide the scripts with the exact settings that we used for this step.

### Small Datasets

For small datasets live BV, use `job_template_simple.sh`,
which just convers and runs epa.

### Larger Datasets

Similar to what we did for the alignment step, we used a parallelization over batches here.
The current version of epa can do better, but it couldn't back when we ran these analyses.
Thus, we just continued to use the batches from running papara.
The setup is the same, so see there for details.

Results
-------------------------

These steps are run for all needed combinations of trees and datasets:

 * `bv_Bact_Constr_backbone`
 * `bv_Bact_Unconstr_backbone`
 * `bv_Gen_Constr_backbone`
 * `bv_Gen_Unconstr_backbone`
 * `hmp_Bact_Constr_backbone`
 * `hmp_Bact_Unconstr_backbone`
 * `hmp_Gen_Constr_backbone`
 * `hmp_Gen_Unconstr_backbone`
 * `tara_Euks_Constr_backbone`
 * `tara_Euks_Unconstr_backbone`
 * `tara_Gen_Constr_backbone`
 * `tara_Gen_Unconstr_backbone`
