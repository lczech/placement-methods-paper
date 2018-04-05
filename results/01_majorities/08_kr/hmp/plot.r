#!/usr/bin/Rscript

library(RColorBrewer)
library(scales)

load("mats.RData")
meta <- read.csv( file="meta_9194_regions.csv", sep=",", head=FALSE )
proj <- read.csv( file="proj.csv", sep=",", head=FALSE )

# Flip the mds components so that the plot gets the same orientation as the Edge PCA one.
# This is valid because the axes and their direction are arbitrary in MDS anyway.
MDS1 <- mds_mat[, 2]
MDS2 <- -mds_mat[, 1]

EdgePC1 <- proj$V2
EdgePC2 <- proj$V3

# Color palettes
col_pal_7 = brewer.pal(7,"Set1")
#col_pal_19 = brewer.pal(19,"Dark2")
#col_pal_7 = rainbow(7)
#col_pal_19 = rainbow(19)
#col_pal_7 = hue_pal()(7)
# col_pal_19 = c(brewer.pal(9,"Set1"), brewer.pal(8,"Dark2"), brewer.pal(3,"Set2")[1])

# col_pal_19[6] = col_pal_19[8]
# col_pal_19[7] = col_pal_19[8]
#
# col_pal_19[9] = col_pal_19[10]
#
# col_pal_19[2] = col_pal_19[3]
# col_pal_19[4] = col_pal_19[3]
# col_pal_19[5] = col_pal_19[3]
#
# col_pal_19[11] = col_pal_19[15]
# col_pal_19[12] = col_pal_19[15]
# col_pal_19[13] = col_pal_19[15]
# col_pal_19[14] = col_pal_19[15]
#
# col_pal_19[16] = col_pal_19[18]
# col_pal_19[17] = col_pal_19[18]

#tmp=col_pal_7[1]
#col_pal_7[1]=col_pal_7[2]
#col_pal_7[2]=tmp

# 0	stool	Stool	GastrointestinalTract
# 1	tongue_dorsum	Mouth (back)	Oral
# 2	palatine_tonsils	Mouth (back)	Oral
# 3	throat	Mouth (back)	Oral
# 4	saliva	Saliva	Oral
# 5	attached_keratinized_gingiva	Mouth (front)	Oral
# 6	hard_palate	Mouth (front)	Oral
# 7	buccal_mucosa	Mouth (front)	Oral
# 8	supragingival_plaque	Plaque	Oral
# 9	subgingival_plaque	Plaque	Oral
# 10	left_retroauricular_crease	Skin	Skin
# 11	right_retroauricular_crease	Skin	Skin
# 12	anterior_nares	Airways	Airways
# 13	left_antecubital_fossa	Skin	Skin
# 14	right_antecubital_fossa	Skin	Skin
# 15	vaginal_introitus	Vagina	UrogenitalTract
# 16	mid_vagina	Vagina	UrogenitalTract
# 17	posterior_fornix	Vagina	UrogenitalTract
# 18	NA	n/a	n/a

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
col_ass_7  <- col_pal_7[as.numeric(cut( meta$V5, breaks = 7))]
col_ass_19 <- col_pal_19[as.numeric(cut( meta$V3, breaks = 19))]

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

# 19 symbols
#########################

# Make PCA graph
svg( "pca_19.svg" )
plot( pca_mat$x[,1:2], pch=symbol, cex=sym_size, col=col_ass_19 )
legend("topright", legend=leg_txt_cat, col=leg_col_cat, pch=leg_sym, cex=leg_size )

# Print MDS
svg( "mds_19.svg" )
plot( MDS1, MDS2, xlim=c(-1.2,1.1), asp = 1, pch=symbol, cex=sym_size, col=col_ass_19 )
legend("topright", legend=leg_txt_cat, col=leg_col_cat, pch=leg_sym, cex=leg_size )

# Print Edge PCA graph
svg( "epca_19.svg" )
plot( EdgePC1, EdgePC2, pch=symbol, cex=sym_size, col=col_ass_19 )
legend("topright", legend=leg_txt_cat, col=leg_col_cat, pch=leg_sym, cex=leg_size )

# 7 symbols
#########################

# Make PCA graph
svg( "pca_7.svg" )
plot( pca_mat$x[,1:2], pch=symbol, cex=sym_size, col=col_ass_7 )
legend("topright", legend=leg_txt_7, col=col_pal_7, pch=leg_sym, cex=leg_size )

# Print MDS
svg( "mds_7.svg" )
plot( MDS1, MDS2, asp = 1, pch=symbol, cex=sym_size, col=col_ass_7 )
legend("topright", legend=leg_txt_7, col=col_pal_7, pch=leg_sym, cex=leg_size )

# Print Edge PCA graph
svg( "epca_7.svg" )
plot( EdgePC1, EdgePC2, pch=symbol, cex=sym_size, col=col_ass_7 )
legend("topright", legend=leg_txt_7, col=col_pal_7, pch=leg_sym, cex=leg_size )
