Overview
-------------------------

The programs in this directory are the prototypes of the clustering methods.

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
