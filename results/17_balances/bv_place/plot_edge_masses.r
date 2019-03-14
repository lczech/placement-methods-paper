#!/usr/bin/Rscript

em <- read.csv( file="edge_masses.csv", sep="\t", head=FALSE, row.names = 1 )

svg( "edge_masses.svg" )
pdf(file="edge_masses.pdf",width=1000,height=1000)
plot(em)
dev.off()
