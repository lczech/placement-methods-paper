Overview
-------------------------

The programs in this directory are the prototypes of the clustering methods.

`bplace_emd_binning`
-------------------------

Testing how binning affects the accuracy and speed of EMD calculations.

`compare_emd_nhd_bplace`
-------------------------

Compare NHD and EMD to each other and to Squash Clustering,
by sorting the distance matrices using the squash cluster merge order.
This yields distance matrices sorted in a way that brings similar rows/colors close to each other,
resulting in a nice visualization of the distnace matrix a as a heat map.

`jplace_emd`
-------------------------

Simple program to calculate the pairwise EMD matrix for a set of jplace files.

`jplace_emd_speed_comp`
-------------------------

Test the speed of our EMD implementation, using different numbers of threads.

`kmeans_jplace`
-------------------------

The program takes four command line arguments:

 1. Number of threads to use.
 2. Number of k clusters to produce.
 3. Input directory with jplace files.
 4. Output directory to write files to.

It then runs kmeans with both masses and imbalances,
and writes the resulting assignments as well as colorized cluster centroid
visualizations.

There is also a version `kmeans_bplace`, which instead of jplace files
expects our internal binary bplace files, for speedup in our testing.

`kmeans_jplace_elbow`
-------------------------

Run kmeans and also output the average distance of the samples to their assigned
centroids, so that we can make elbow plots.
