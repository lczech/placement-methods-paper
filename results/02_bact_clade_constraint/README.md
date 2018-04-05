Overview
-------------------------

Test how the accuracy of the ART changes 
when we do a high level tax constraint to get our five test clades constrained.
This is conducted for the Bacteria tree, and uses the data from `11_consensus_seqs`.

To this wend, we trimemd the taxonomy to just contain our five clades of interest,
see `00_reference` for the scripts and resulting alignment.
The Silva sequences of the five clades were then placed as usual,
and finally evaluated in terms of "how many of them went into the correct clade"
in `03_eval`, which uses the program `code/art/silva_subclade_hits.cpp`.
