#!/usr/bin/Rscript

library( MASS )
library( cluster )

input_name <- "emd"
input_file <- paste(input_name, ".mat", sep = "")
dist_data  <- read.table( input_file, sep=" " )
dist_mat   <- do.call( rbind, dist_data )

clusters <- pam( dist_mat, 5, diss=TRUE, keep.diss=TRUE )

svg( paste("clusplot_", input_name, ".svg" , sep = ""))
# clusplot( clusters )
# clusplot( clusters, lines=0, col.p=c("dark green", "blue", "red"))
# clusplot( clusters, lines=0, col.p=clusters$cluster )
# clusplot( clusters, lines=0, col.p=c("red", "green", "blue")[clusters$cluster] )

# palette( heat.colors(6) )


meta <- read.csv( file="meta_all", sep="\t", head=FALSE )
# col_fac <- factor( meta$V2 )

rbPal <- colorRampPalette(c('blue','red'))
colpal <- rbPal(6)[as.numeric(cut( meta$V2, breaks = 6))]
cuts <- c( "0", "2", "4", "6", "8", "10" )

# par(pch=19)

# summary(col_fac)
# clusplot( clusters, lines=0, col=col_fac )
# clusplot( clusters, col.p=col_fac )
clusplot( clusters, col.p=colpal, plotchar=FALSE, col.pch=0  )
# clusplot( clusters )

legend("topright", cuts, col=rbPal(6), pch=19 )

# colorthings=scan("col")
# # colorthings=1
# clusplot( clusters, lines=0, col.p=colorthings )

# clusplot( clusters, lines=0, col.p=clusters$cluster )
dev.off()



fit <- isoMDS( dist_mat, k=2 )
# fit <- cmdscale(dist_data, eig = TRUE, k = 2)

x <- fit$points[, 1]
y <- fit$points[, 2]

svg( paste( "isomds_", input_name, ".svg" , sep = ""))
plot( x, y, col=colpal, pch=19 )
legend("topright", cuts, col=rbPal(6), pch=19 )

# plot(x, y)
# text(x, y)
dev.off()

# write( colpal, "map_emd_5.cluster", ncolumns=1 )

# library(igraph)
# g <- graph.full(nrow(dist_data))
# layout <- layout.mds(g, dist = as.matrix(dist_data))
# plot(g, layout = layout, vertex.size = 3)
