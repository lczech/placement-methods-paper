Overview
-------------------------

The programs in this directory were used for preparing the raw sequence data,
so that the methods and other programs can work with them,
as well as analyzing certain aspects of the data.
The programs use genesis v0.19.0

`hmp_histograms.cpp`
-------------------------

Build histograms of sequence lengths for the HMP data.
This was used to estimate which sequences we want to consider too short or too long
for our anlysis, that is, what we want to consider a "typical" sequence for this dataset.
See `data/hmp/README.md` for the results.

`silva_tree_eval_blacklist.cpp`
-------------------------

Write a blacklist file that lists all sequence names (SEQ_xxxxxx_) 
that we want to try to exclude from evaluation.
This evaluation is mentioned in the ART supplement, where we wanted to assess
the effect of excluding "bad" sequences. We consider as bad:

 * all sequences form the Sativa mislabel list
 * all sequences that contain "incertae", "unclassified" or "unknown" in their taxopath

This list has 25,910 entries out of 598,470 sequences, or 4.3% of the data.
The results can be found in `data/silva/blacklist.txt`.
