Overview
-------------------------

Here, we provide the experiments on the BV data using the reference sequences
from the original publication

> "Bacterial communities in women with bacterial vaginosis: High resolution phylogenetic analyses reveal relationships of microbiota to clinical criteria".
> Srinivasan et al., PLoS One, vol. 7, no. 6, p. e37818, Jan. 2012.
> http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037818

We used these data for testing our methods against the established analysis methods.

 * `00_data`: The reference data (sequences and tree) as provided in the online supplement "Reference Package S1".
   Unfortunately, we cannot also share their query sequences here, as those are patient data which we
   obtained directly from Sujatha Srinivasan
 * `01_treesearch`: We re-inferred the reference tree from their alignment, just to make sure
   that we use the same model parameter for the tree that are later used for the placement.
   All the following steps were thus conducted on both, their original and our re-inferred tree.
   For the final paper, we used our tree.
 * `02_align`: Aligning the queries against the reference.
 * `03_epa`: Placing the sequences on the tree.
 * `04_dists`: Some tests of measuring and visualizing distances between the samples.
 * `05_emd`: Plots of Edge PCA, MDS and PCA on the dataset, as shown in our paper.
 * `05_epca`: Running Edge PCA on the data, thereby re-creating a figure from the original paper.
 * `05_epca_new`: Final run of Edge PCA using our new implementation in genesis/gappa, which
   uses the color scheme that we used throughout our paper.
 * `06_kmeans`: Test runs of kmeans on the dataset.
 * `08_squash_kmeans`: Comparison of Squash Clustering to kmeans. Here, we create the plots 
   of squash cluster trees colored with the kmeans cluster assignments.
 * `09_squash_nugent`: Coloring of the squash cluster tree with the nugent score per sample.
