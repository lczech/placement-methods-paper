Overview
-------------------------

This directory contains the scripts of aligning the query sequences 
to the reference consensus sequences and the reference tree inferred from them.
We use PaPaRa for this step.

We do not offer all alignments here online, as this would be too much data.
Instead, we provide the scripts with the exact settings that we used for this step.

### Small Datasets

For small datasets live BV, use `job_template_simple.sh`,
which just runs PaPaRa.

### Larger Datasets

Currently, PaPaRa is not well paralelized and thus a bottleneck in the analyses.
Hence, we created our own parallel version of it, so that at least we can use
differnt MPI nodes at once.
This version can be found here: https://github.com/lczech/papara_nt
In order to use it, the data has to be split into chucks, which is already done
by our sequence preprocessing.

Here, we now need to make sure to run the chunks correctly and provide papara
with a sufficient number of input files.

To achive this, frist we create a list of the chunks that are availabe,
and split this list into batches of, say, 20 chunks, which reflects our MPI
cluster. That is, each MPI run of our papara version gets 20 chunks to run
on 20 nodes. To create this, navigate to the chunk directory of the sequence dataset
and run:

    ls chunks/ > batches.txt
    mkdir batches
    split -l 20 batches.txt batches/batch_
    
This data is then used by the `create_jobs.sh` script to create a job with 20 MPI nodes
for each batch. The script creates copies of `job_template.sh`, and replaces
the necessary placeholders for the batches and chunks.

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
