Overview
-------------------------

This directory contains the code and scripts that we used for the papers.
That is, using this, one should be able to reproduce our exact results
as presented in the papers.
We modified our original code by removing hard coded file paths and
similar details that would have made reproduction of results harder.
By this, we tried to make the code as portable as feasible.
However, it is still possible that some passages have to be adjusted
for other environments in order to work properly.

All programs here are in C++ and intended to be used with our library
[genesis v0.19.0](https://github.com/lczech/genesis/releases/tag/v0.19.0).
The easiest way to set the programs is to use 
the [Apps feature](http://doc.genesis-lib.org/setup.html#setup_apps)
of genesis, where you simply copy the `cpp` files from here into the `apps`
directory of genesis and compile them.

 * `art`: Our prototype implementation of the Automatic Reference Tree method,
   as well as the code that we used for the testing and evaluation.
 * `clustering`: Prototypes of kmeans clustering, and our comparsion
   with Squash Clustering.
 * `data`: Data preprocessing programs that are not specific to one method,
   but were used to get the sequence data into usable formats in the first place.
 * `multilevel`: Prototypes for the multilevel placement approach,
   that is, chunkify and unchunkify, extract, etc.
 * `tests`: Some old code and scripts that we used during the development process.
   These are just provided as-is, and probably won't compile without some changes.
   We keep them here just in case we later need to revisit some of the
   early details of our methods.
 * `tools`: Some useful programs that might be reusable beyond our papers as well,
   for example merging jplace files or converting sequence file formats.
 * `visualization`: Prototypes of the visualization methods, e.g, 
   Edge Dispersion and Edge Correlation,as well as related programs.
