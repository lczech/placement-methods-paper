Overview
-------------------------

The programs in this directory were used for preparing the raw sequence data,
so that the methods and other programs can work with them,
as well as analyzing certain aspects of the data.
The programs use genesis v0.19.0

`fasta_cleanup`
-------------------------

Simple tool to prepare the silva sequences. Replace chars: "." --> "-"

`hmp_histograms`
-------------------------

Build histograms of sequence lengths for the HMP data.
This was used to estimate which sequences we want to consider too short or too long
for our anlysis, that is, what we want to consider a "typical" sequence for this dataset.
See `data/hmp/README.md` for the results.

`silva_fasta_filter`
-------------------------

Tool for extracting a part of the Silva sequences belonging to one taxonomic path.

