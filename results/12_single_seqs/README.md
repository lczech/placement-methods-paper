Overview
-------------------------

This is the test for Silva single sequence representation of our entropy pruned trees.

We take the same entropy pruned taxa as we did for the main part of the paper (backbone tree),
but instead of representing each taxon by the consensus sequence of all sequences that
belong to that taxon, we represent it by only one sequence.
This sequence is the closest of all sequences to the consensus (in terms of different nucleotides), 
in order to also be a good representative of that taxon. But it does not contain ambiguity chars,
and also might have gaps in different places than the consensus sequence.

The idea is that this single sequence does not represent its clade as well as
a consensus sequence, because it does not incorporate as much biological information.
However, it might also happen that by using consensus, we are blurring the signal,
and thus get worse results.
This directory contains the tests to determine which of those is the case.
