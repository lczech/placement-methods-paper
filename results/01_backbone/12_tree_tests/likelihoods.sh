#!/bin/bash

for r in "Arch" "Bact" "Euks" "Gen"; do
    for t in "Constr" "Unconstr"; do
        echo "Eval: ${r} ${t}"
        raxml -f n -z 01_trees/${r}_${t}_backbone/best_tree.newick -s 00_reference/${r}/Ref.phylip -m GTRGAMMAX -n ${r}_${t} -T 4 &> ${r}_${t}.log
        grep "Tree 0 Likelihood" ${r}_${t}.log
        echo
    done
done
