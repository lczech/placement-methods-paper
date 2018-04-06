Overview
-------------------------

The programs in this directory are the prototypes of the visualization methods.

`dispersion_trees`
-------------------------

The program takes two command line arguments:

 1. Input directory with jplace files.
 2. Output directory to write files to.

It then calculates several types of dispersion measures and visualizes
them on the tree, writing one colored tree per measure, both in svg and nexus format.

`jplace_squash`
-------------------------

Run Squash Clustering on a set of jplace files:

 1. Number of threads to use.
 2. Input directory with jplace files.
 3. Output directory to write files to.

The program procudes a squash cluster tree in newick format.

`jplace_squash_nugent`
-------------------------

Special version of the normal squash clustering that additionally writes
an svg file with the cluster tree where the tip nodes are annotated
with circles in colors that represent the nugent score of their samples.

`label_matrix`
-------------------------

Turn the cluster assignment of kmeans and the original labelling
into a heatmap showing which clusters contain how many samples from each label.

