#!/usr/bin/Rscript

library(RColorBrewer)
library(scales)

# Calc mds and pca
# dists <- read.csv( file="balance_distances.csv", sep="\t", head=FALSE )
# pca_mat <- prcomp(dists)
# mds_mat <- cmdscale(dists)
# save( pca_mat, mds_mat, file = "balance_distances.RData" )


proj <- read.csv( file="factor_balances_all.csv", sep=",", head=FALSE )

X <- proj$V2
Y <- proj$V3

# Color palettes
col_pal_7 = brewer.pal(7,"Set1")

# Summarize certain body sites, in order to get a clear plot.
# The order of entries in this list is important, as it is the order of the category labels 0-19.
# We hand select the order of colors, so that nearby clusters get different colors,
# in order to be distinguishable from each other.
# Also, we select them to be somehow relatable to the body site.
set = brewer.pal(9,"Set1")
col_pal_19 = c(
    set[7],	#  1 Stool
    set[5],	#  2 Mouth (back)
    set[5],	#  3 Mouth (back)
    set[5],	#  4 Mouth (back)
    set[1],	#  5 Saliva
    set[4],	#  6 Mouth (front)
    set[4],	#  7 Mouth (front)
    set[4],	#  8 Mouth (front)
    set[3],	#  9 Plaque
    set[3],	# 10 Plaque
    set[2],	# 11 Skin
    set[2],	# 12 Skin
    set[9],	# 13 Airways
    set[2],	# 14 Skin
    set[2],	# 15 Skin
    set[8],	# 16 Vagina
    set[8],	# 17 Vagina
    set[8],	# 18 Vagina
    0 # set[6]	# 19 n/a
)

# Set1 Colors
# 1 red
# 2 blue
# 3 green
# 4 purple
# 5 orange
# 6 yellow
# 7 brown
# 8 pink
# 9 gray

# Make a legend. The order here is arbitrary (we order head to toe, kind of...),
# but needs to be maintainted between the two lists.
leg_txt_cat = c(
    "Mouth (back)",
    "Mouth (front)",
    "Saliva",
    "Plaque",
    "Airways",
    "Skin",
    "Stool",
    "Vagina"
    # "n/a"
)
leg_col_cat = c(
    col_pal_19[2],
    col_pal_19[6],
    col_pal_19[5],
    col_pal_19[9],
    col_pal_19[13],
    col_pal_19[11],
    col_pal_19[1],
    col_pal_19[16]
    # col_pal_19[19]
)

# For Set1, we flip yellow and brown, so that the interesting data is better visible
tmp=col_pal_7[6]
col_pal_7[6]=col_pal_7[7]
col_pal_7[7]=tmp

# Color assignments
col_ass_7  <- col_pal_7[as.numeric(cut( proj$V19, breaks = 7))]
col_ass_19 <- col_pal_19[as.numeric(cut( proj$V18, breaks = 19))]

# http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
symbol <- 20
sym_size <- 0.4
leg_sym <- 19
leg_size <- 1

# Legend
#########################

#stool	0
#mouth	1
#ear	2
#nose	3
#arm	4
#vagina	5
#NA	6

leg_txt_7=c("Stool","Mouth","Ear","Nose","Arm","Vagina","n/a")
# leg_txt_19=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18")

# 7 symbols
#########################

# Print MDS
svg( "pf_all.svg" )
#plot( X, Y, asp = 1, pch=symbol, cex=sym_size, col=col_ass_7 )
#legend("topright", legend=leg_txt_7, col=col_pal_7, pch=leg_sym, cex=leg_size )


plot( X, Y, pch=symbol, cex=sym_size, col=col_ass_19 )
legend("topright", legend=leg_txt_cat, col=leg_col_cat, pch=leg_sym, cex=leg_size )


