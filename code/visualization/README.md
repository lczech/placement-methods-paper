Overview
-------------------------

The programs in this directory are the prototypes of the visualization methods.

`cluster_tree_metadata`
-------------------------

Take a cluster tree (from Squash clustering), as well as a metadata file, and annotate the
leaves of the tree with colors according to metadta feature values.
Because each leaf of the cluster tree corresponds to a samples,
there is one leaf per metadata table row. The result then visualizes the metadata distribution
of samples on the clsuter tree, as for example shown in the ART evaluation of the BV data,
using their Nugent score metadata.

`correlation_trees`
-------------------------

Prototype implementation of the Edge Correlation method.
Takes a directory with jplace files, a metadata table in csv format, and an ouput dir.

`dispersion_trees`
-------------------------

The program takes two command line arguments:

 1. Input directory with jplace files.
 2. Output directory to write files to.

It then calculates several types of dispersion measures and visualizes
them on the tree, writing one colored tree per measure, both in svg and nexus format.

`jplace_epca_vis`
-------------------------

Re-implementation and our prototype of the Edge PCA visualization.
Instead of using branch width to visualize the eigenvalues, we use colors.
This program was used for the Edge PCA visualizations in the paper.

`jplace_squash`
-------------------------

Run Squash Clustering on a set of jplace files:

 1. Number of threads to use.
 2. Input directory with jplace files.
 3. Output directory to write files to.

The program procudes a squash cluster tree in newick format.

`jplace_squash_kmeans`
-------------------------

Visualization of kmeans cluster assignments as circles 
on the tips of the Squash Cluster tree.
Use for the tests and the repsective figure of the paper.

`jplace_squash_nugent`
-------------------------

Special version of the normal squash clustering that additionally writes
an svg file with the cluster tree where the tip nodes are annotated
with circles in colors that represent the nugent score of their samples.

`label_matrix`
-------------------------

Turn the cluster assignment of kmeans and the original labelling
into a heatmap showing which clusters contain how many samples from each label.
Used for the heatmap figures in the paper.
