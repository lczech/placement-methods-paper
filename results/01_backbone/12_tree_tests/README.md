Likelihoods
-------------------------

Provided in`likelihoods.sh`: Evaluate the likelihood scores of the automatic reference trees, 
particularly of the constrained ones.
Result: the unconstr have a better likelihood than the constrained ones.

RF Distances
-------------------------

Mini test to get the rf distances between the unconstrained and constrained trees.

Combine both trees into one file:

    $ cat Arch_*/best_tree.newick > Arch_trees.newick
    $ cat Bact_*/best_tree.newick > Bact_trees.newick
    $ cat Euks_*/best_tree.newick > Euks_trees.newick
    $ cat Gen_*/best_tree.newick > Gen_trees.newick

Calculate distances:

    $ raxml -m GTRCAT -f r -z Arch_trees.newick -n Arch_trees > Arch.log
    $ raxml -m GTRCAT -f r -z Euks_trees.newick -n Euks_trees > Euks.log
    $ raxml -m GTRCAT -f r -z Bact_trees.newick -n Bact_trees > Bact.log
    $ raxml -m GTRCAT -f r -z Gen_trees.newick -n Gen_trees > Gen.log

The svg files are created with `code/tools/newick_compare.cpp` and show differing splits.
This was used to state that the differences are mainly on the inner branches.

Significance Tests
-------------------------

We used IQ-tree to run significance tests comparing constrained vs unconstrained trees.
The script is `significance.sh`, results are shown in `significance.txt`.
They use the combined trees from the RF distnaces, see above,
as well as the alignments from `00_reference`.
