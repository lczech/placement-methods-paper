placement-methods-paper
=========================

Accompanying code repository with scripts, programs and data for our papers:

> Methods for Inference of Automatic Reference Phylogenies and Multilevel Phylogenetic Placement.<br />
> Lucas Czech, Pierre Barbera, and Alexandros Stamatakis.<br />
> Bioinformatics, 2018. https://doi.org/10.1093/bioinformatics/bty767<br />
> <!-- bioRxiv, 2018. https://doi.org/10.1101/299792 -->

and

> Scalable Methods for Analyzing and Visualizing Phylogenetic Placement of Metagenomic Samples.<br />
> Lucas Czech and Alexandros Stamatakis.<br />
> PLOS One, 2019. https://doi.org/10.1371/journal.pone.0217050<br />
> <!-- bioRxiv, 2019. https://doi.org/10.1101/346353 -->

Structure
-------------------------

The repository is structed as follows:

 * `code`: C++ code to be used with our [genesis](https://github.com/lczech/genesis) library. 
   The directory contains all prototypes and helper programs that were used to
   evaluate the data and make the plots reported in the papers.
 * `data`: Contains information and some processed files about the empirical 
   datasets that were used as a basis for the analyses and evaluations.
 * `results`: A large collection of scripts, data, results and figures that were used in the papers.
   The directory contains all analyses that we conducted for the papers.
 * `software`: Information about the external software that was used.

The purpose of this repository is to allow reconstruction of the analyses that we ran for the papers.
If you are however interested in running analyses on your own data,
that is, use our methods, have a look at [gappa](https://github.com/lczech/gappa).

**IMPORTANT REMARK:** During the submission of the paper on Automatic Reference Phylogenies,
we decided to rename the method from ART to PhAT.
While this renaming is now consistently applied to both the paper and the gappa implementation,
in this repository, we did not rename it. Thus, be aware that all references to ART
are meant to refer to PhAT instead.

Quick Links
-------------------------

The file [`results/msas_and_trees.zip`](results/msas_and_trees.zip)
is a shortcut for users who are simply interested in the alignments and trees produced 
by our PhAT/ART method for their own analyses.
The file contains MSAs, taxonomies, constrained and unconstrained trees
for all four domains that we used in the paper (General, Archaea, Bacteria, Eukaryota).
The files are identical to the ones found in [`results/01_backbone`](results/01_backbone),
but renamed for simplicity.
