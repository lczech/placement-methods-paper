This is a short test to see how many of the placed sequences belong to a collection of "hard to place" genera, 
according to http://jcm.asm.org/content/45/9/2761.short

For this, we use the table at http://jcm.asm.org/content/45/9/2761/T2.expansion.html
and scan the genera and species listed there in our data, using the script here.

then, we use the two resulting counts files, one for genera one for species, 
to get to total number of hard to classify sequences that are part of those genera and species.
that is:

3264 species sequences
95553 genera sequences


finally, we filter these results again to find the number of corresponding taxa in the reference tree 
that belong to those genera, using find_mean_taxa

for this, we copied the tips taxa lists from 02_multilevel/00_reference/subtree_aln
and scan them for entries as well

224 taxa_filtered
are the taxa of the mean genera


Subtree Phyla
==========================

Cyanobacteria
94 taxa
11,379 sequences

Proteobacteria
1,569 taxa
200,083 sequences

Firmicutes
652 taxa
138,636 sequences

Bacteroidetes
440 taxa
49,214 sequences

Actinobacteria
494 taxa
51,160 sequence

total: 450,472

^^^ not correct, this is all taxa including higher level ones.
if we just count tips, we get

   406 Bacteria_Actinobacteria.tips.tax
   406 Bacteria_Bacteroidetes.tips.tax
    79 Bacteria_Cyanobacteria.tips.tax
   583 Bacteria_Firmicutes.tips.tax
  1372 Bacteria_Proteobacteria.tips.tax
  2846 total

which is also the count of the fasta files used for inferring the trees, so that looks more correct!

of these, 224 are mean ones according to the paper
of the sequences,

1572    genus_details/Aeromonas_filtered
12542   genus_details/Bacillus_filtered
389 genus_details/Bordetella_filtered
25837   genus_details/Burkholderia_filtered
4267    genus_details/Campylobacter_filtered
147 genus_details/Edwardsiella_filtered
24431   genus_details/Enterobacter_filtered
2263    genus_details/Neisseria_filtered
15210   genus_details/Pseudomonas_filtered
8895    genus_details/Streptococcus_filtered
95553

are in mean genera

number of sequences:
02_sequences/subtree_sequences/subtree_sequences_no_silva_gaps>  grep ">" * | wc -l

yields: 450313

