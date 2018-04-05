#!/usr/bin/Rscript

# Read pairwise emd dists mat and edge pca projection
kr_dists <- read.csv( file="emd.mat", sep=" ", head=FALSE )
meta <- read.csv( file="meta_9194.csv", sep=",", head=FALSE )

# Run PCA and MDS
pca_mat <- prcomp(kr_dists)
mds_mat <- cmdscale(kr_dists)

# Use three colors for three clusters
#rgb_pal <- colorRampPalette(rainbow(19))
pal <- rainbow(19)[as.numeric(cut( meta$V3, breaks = 19))]

# Make PCA graph
svg( "pca.svg" )
plot( pca_mat$x[,1:2], pch=19, col=pal )

# Print MDS
MDS1 <- mds_mat[, 1]
MDS2 <- mds_mat[, 2]

svg( "mds.svg" )
plot( MDS1, MDS2, asp = 1, pch=19, col=pal )

# Print Edge PCA graph
EdgePC1 <- -proj$V2
EdgePC2 <- -proj$V3

svg( "epca.svg" )
plot( EdgePC1, EdgePC2, pch=19, col=pal )

