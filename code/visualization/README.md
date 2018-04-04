Overview
-------------------------

The programs in this directory are the prototypes of the visualization methods.
The programs use genesis v0.19.0

`dispersion_trees.cpp`
-------------------------

The program takes two command line arguments:

 1. Input directory with jplace files.
 2. Output directory to write files to.

It then calculates several types of dispersion measures and visualizes
them on the tree, writing one colored tree per measure, both in svg and nexus format.
