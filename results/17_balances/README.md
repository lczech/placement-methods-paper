Overview
-------------------------

The directory contains scripts and result files of our adaptations
of the PhyILR transformation and Phylofactization to phylogenetic placement data.
Many of the subdirectories contain dirs named `tw` and `no_tw`, which stand
for the two variants with taxon weighting and without (no) taxon weighting.

 * `bv_place`: Evaluation of placeing the BV dataset on their ref tree, then
   running our adapations of PhILR and Phylofactorization on it.
 * `bv_swarm_consensus`: Evaluation of the BV dataset with the original Phylofactorization,
   when using swarm for OTU clustering.
 * `bv_vserach_cluster`: Same, but for clustering with vsearch instead.
 * `hmp_pf_all`: Eval of the full HMP dataset with our adapations.
 * `hmp_pf_of_600`: Only using the 600 oral and fecal samples of the hmp dataset instead.
