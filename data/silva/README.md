Overview
-------------------------

Silva sequence database.

> "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools".
> Quaste et al., Nucleic Acids Res., vol. 41, no. D1, pp. D590–D596, Jan. 2013.
> https://www.arb-silva.de/

> "The SILVA and ‘All-species Living Tree Project (LTP)’ taxonomic frameworks".
> Yilmaz et al., Nucleic Acids Res., vol. 42, no. D1, pp. D643–D648, Jan. 2014.

The data was downloaded from https://www.arb-silva.de/no_cache/download/archive/release_123_1/

 * `SILVA_123.1_SSURef_Nr99_tax_silva_full_align_trunc.fasta.gz`: 598 470 aligend sequences. 
   This is the main data used for our paper, that is, these sequences were used to generate
   consensus sequences.
 * `tax_slv_ssu_123.1.txt`: 11 860 taxonomic ranks of the database.
   This is the taxonomy that we used for the ART method and its entropy calculations.
   
The data was processed with the programs in `code/art`. See there for more on this.

Additional Files
-------------------------

### `mislabels.txt`

Extracted from the file `nar-00648-z-2016-File011.xlsx` in the Sativa supplement file at
[https://academic.oup.com/nar/article/44/11/5022/2468320/Phylogeny-aware-identification-and-correction-of](https://academic.oup.com/nar/article/44/11/5022/2468320/Phylogeny-aware-identification-and-correction-of).
It contains the sequence names of the sequences that were considered mislabels by Sativa.
It is used as input to the program `code/data/silva_tree_eval_blacklist.cpp`.

### `blacklist.txt`

This file was created using the mislabels as well as certain keywords that signal
that a sequence is not trustworthy.
We used the Genesis program `code/data/silva_tree_eval_blacklist.cpp` for this.
It lists all sequences with "incertae" "unclassified" or "unknown" in their name,
and all mislabel sequences.
The format is: Our name SEQ_xxxxxx, followed by the reason it appears in the blacklist
(either one of the three keywords, or the accession name as used in the mislabel list).

### `silva_full_align_reduced.gaps`

This file contains bit-wise information about which sites in the full silva alignment are all gaps.
It was created by the program in `code/tests/silva_gap_finder.cpp`, 
and can be use by `code/tests/silva_gap_user.cpp` in order to remove those gap sites from 
other files based on the silva alignment.
We used this to make some of our files smaller.
