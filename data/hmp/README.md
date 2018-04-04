Overview
-------------------------

Human Microbiome Project dataset.

> "Structure, function and diversity of the healthy human microbiome,".
> Huttenhower et al., Nature, vol. 486, no. 7402, pp. 207–214, Jun. 2012.

> "A framework for human microbiome research".
> Methé, et al.,  Nature, vol. 486, no. 7402, pp. 215–221, Jun. 2012.

The data was obtained from their Phase 1 ([https://www.hmpdacc.org/hmp/](https://www.hmpdacc.org/hmp/)) in June 2017.
We used the `HM16STR` sequences from [http://hmpdacc.org/HM16STR/all/](http://hmpdacc.org/HM16STR/all/),
which we downloaded from `public-ftp.hmpdacc.org/HM16STR`.

The processing of the datasets done by the HMP project is described at 
[ftp://public-ftp.hmpdacc.org/HM16STR/README.txt](ftp://public-ftp.hmpdacc.org/HM16STR/README.txt).

We used two sub-datasets:

## Dataset `published`

We stored them in directory `by_sample`.
It contains the published reads, all healthy,
download links provided in `HM16STR.csv`, downloaded from http://hmpdacc.org/HM16STR/

## Dataset `healthy`

We stored them in directory `SRP002860`.
It contains additional reads, not published in studies, all healthy,
download links provided in `HM16STR_healthy.csv`, downloaded from http://hmpdacc.org/HM16STR/healthy/

Processing
-------------------------

First, we used the `make_links.sh` script to obtain a directory `links`,
which contains links to all (except a few overlength) sequences of
`by_sample` and `SRP002860`.

There are this many sequences per primer in total:

    count V1-V3  41,498,877
    count V3-V5  76,180,530
    count V6-V9   1,023,560
    count all   118,702,967

Then, we use the program in `code/data/hmp_histograms.cpp` to obtain length histograms for those
sequences, splitted by primer, in order to decide which min max length
we want to use for filtering.

These histograms are prepared in `length_histograms.ods` for easy 
overview. The resulting data is used to decide on the min and max
boundaries for valid sequence length.

We decided to use a min of 150, and following max lengths:

    primer=V1-V3 --> max 540
    primer=V3-V5 --> max 575
    primer=V6-V9 --> max 540

We set the max to the position where the last declining number appears.
That is, the next count in the histogram is higher than the previous one.

Lastly, we used the program `code/art/filter_chunkify.cpp` to filter and split into chunks,
which will create a directory `filtered` for the filtered fasta files,
and a directory `chunks` for chunks of 50,000 sequences each, 
and a `maps` directory for the abundance map files.
The program is a prototype of the `gappa chunkify` command.

This results in a total of 104,941,453 sequences in 2,099 chunks.

Those chunks are then the input into our pipeline. Each sequence
label in there has its original fasta file name as prefix, so that we
later can split it again.

Additional Files
-------------------------

`categories.csv`, `label_list.txt`, `label_body_site.txt`, `label_body_site_list.txt` and `meta.csv`
are representations of the body site label meta-data feature of the dataset.

`meta_9194.csv` contains the sorted meta entries for the filtered list of samples 
that we actually used for our test. created by:

    $ ls filtered/ | sed "s/\.fsa//g" | sort | xargs -I LINE grep LINE meta/meta.csv | sed "s/\t/,/g" > meta/meta_9194.csv
    

`meta_9194_overview.ods` gives an overview of how many samples come from each body site.

`meta_9194_regions.csv`, `label_body_site_list_8.txt` and `label_body_site_list_18.txt`
were used for evaluating the k-means clustering using different values of k.

Tools
-------------------------

Get how many samples there are for each body site:

    cat meta_9194_regions.csv | cut -d, -f 2 | sort | uniq -c
    
    541 anterior_nares
    600 attached_keratinized_gingiva
    597 buccal_mucosa
    566 hard_palate
    290 left_antecubital_fossa
    596 left_retroauricular_crease
    298 mid_vagina
      2 NA
    599 palatine_tonsils
    301 posterior_fornix
    328 right_antecubital_fossa
    604 right_retroauricular_crease
    529 saliva
    600 stool
    595 subgingival_plaque
    608 supragingival_plaque
    638 throat
    610 tongue_dorsum
    292 vaginal_introitus
    
Might be useful.
