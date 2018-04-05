Mini test to get the rf distances between trees of difference consensus methods.

Combine the respective trees into one file:
cat 11_consensus_seqs/01_trees/trees/*_Archaea_Unconstr_best_tree.newick > Archaea_trees.newick
cat 11_consensus_seqs/01_trees/trees/*_Bacteria_Unconstr_best_tree.newick > Bacteria_trees.newick
cat 11_consensus_seqs/01_trees/trees/*_Eukaryota_Unconstr_best_tree.newick > Eukaryota_trees.newick
cat 11_consensus_seqs/01_trees/trees/*_General_Unconstr_best_tree.newick > General_trees.newick

Calculate distances:
raxml -m GTRCAT -f r -z Archaea_trees.newick -n Archaea_trees > Archaea.log
raxml -m GTRCAT -f r -z Bacteria_trees.newick -n Bacteria_trees > Bacteria.log
raxml -m GTRCAT -f r -z Eukaryota_trees.newick -n Eukaryota_trees > Eukaryota.log
raxml -m GTRCAT -f r -z General_trees.newick -n General_trees > General.log

Get all distances in one file:
cat RAxML_RF-Distances.* > all_dists

this gives RF distances across all trees and consensus methods:
min 0.29037
max 0.738189
avg 0.495003125

