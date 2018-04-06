Overview
-------------------------

The programs in this directory are the prototypes of the multilevel methods.

`dispersion_trees`
-------------------------

The program takes two command line arguments:

 1. Abundance map file directory
 2. Jplace directory with the chunified jplace files.
 3. Output directory to write files to.

It then undoes the chunkfiy process, that is, it writes per-sample jplace files
with their abundances per sequence.
