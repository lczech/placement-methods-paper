#!/usr/bin/Rscript

proj <- read.csv( file="factor_balances.csv", sep="\t", head=FALSE )

# head(proj)

# palette( heat.colors(6) )
# palette( topo.colors(6) )
# palette( terrain.colors(6) )

rbPal <- colorRampPalette(c('blue','red'))
colpal <- rbPal(6)[as.numeric(cut( proj$V4, breaks = 6))]

# cuts<-levels(cut( proj$V4,breaks = 6 ))
# cuts<-gsub(","," - ",cuts)
# cuts<-gsub("\\(","[",cuts)
cuts <- c( "0", "2", "4", "6", "8", "10" )

svg( "phylo_factors.svg" )
# plot( proj$V2, proj$V3, pch=19, col=factor( proj$V4  ))
plot( proj$V2, proj$V3, pch=19, col=colpal )

legend("topright", cuts, col=rbPal(6), pch=19 )
