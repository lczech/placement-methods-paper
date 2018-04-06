Overview
-------------------------

The directory contains the full analysis pipeline of the backbone trees for our ART method:

 * `00_reference`: The raw results of the method, that is, consensus sequences for taxonomic ranks,
   built from the Silva sequences and taxonomy.
 * `01_trees`: The sequences and the original taxonomy were used to infer constrained and
   unconstrained trees.
 * `02_sequences`: Shortcut to the query sequences to be placed on those trees for eval purposes.
 * `03_align`: Align the query sequences of our test datasets to the consensus sequence alignments and trees.
 * `04_epa`: Run the placement algorithm to get jplace files.
 * `06_samples`: As we chunkified the data before placing it, we need to unchunkify in order to get per-sample
   placement files.
 * `08_emd`: Some test runs of calculating KR distance (aka EMD) as well as first
   evaluations of Edge Dispersion.
 * `09_viz`: Some tests that show the per-edge masses visualized on the trees.
 * `09_viz_correlation_trees`: Visualizations of Edge Correlation for the paper, using BV and Tara data.
 * `09_viz_dispersion`: Visualizations of Edge Dispersion of the BV data for the paper.
 * `10_kmeans`: k-means runs and visualizations on the BV and HMP datasets.
 * `11_tree_eval`: Evaluation of the ART trees, by placing the Silva sequences on them and measuring
   distances to expected placement branches.
 * `12_tree_tests`: Likelihood and significance tests, RF distance etc of the ART trees.
