#!/usr/bin/Rscript

# Read pairwise emd dists mat and edge pca projection
kr_dists <- read.csv( file="emd.mat", sep=" ", head=FALSE )

# Run PCA and MDS
pca_mat <- prcomp(kr_dists)
mds_mat <- cmdscale(kr_dists)

save( pca_mat, mds_mat, file = "mats.RData" )
