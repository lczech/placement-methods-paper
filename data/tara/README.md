Overview
-------------------------

Tara Oceans dataset.

> "A holistic approach to marine Eco-systems biology".
> Karsenti et al., PLoS Biol., vol. 9, no. 10, pp. 7–11, 2011.

> "Structure and function of the global ocean microbiome".
> Sunagawa et al., Science (80-. )., vol. 348, no. 6237, pp. 1–10, 2015.

> "Plankton networks driving carbon export in the oligotrophic ocean".
> Guidi et al., Nature, vol. 532, no. 7600, p. 465—470, Apr. 2016.

We downloaded the data from [https://www.ebi.ac.uk/ena/data/view/PRJEB6610](https://www.ebi.ac.uk/ena/data/view/PRJEB6610).
See there for a description of the dataset.
Further information can also be found at [http://ocean-microbiome.embl.de/companion.html](http://ocean-microbiome.embl.de/companion.html).

Download
-------------------------

Ftp links for downloading as well as the scripts 
that we used for automating the download process
are found in directory `download`.

When executed, the scripts create a directory `fastq` and 
fills it with the raw fastq sequence data from the downloads.
These files then were processed as follows.

Data Processing
-------------------------

Scripts and tools for processing the raw fastq data are provided in directory `processing`.
The processing needs three programs: `pear`, `cutadapt`, and `vsearch`. 
See `software` for links and information about them.

First, the fastq files are assembled via `assemble_paired_reads.sh`.
Then, some cleaning steps are performed in `primer_clipping_quality_extraction_dereplication.sh`,
which generates fasta files for the samples.
Finally, `recompress.sh` helps to save some storage.

Lastly, we used the program `code/art/filter_chunkify.cpp` to filter and split into chunks,
which will create a directory `filtered` for the filtered fasta files,
and a directory `chunks` for chunks of 50,000 sequences each, 
and a `maps` directory for the abundance map files.
The program is a prototype of the `gappa chunkify` command.

Meta-Data
-------------------------

Meta-data are provided in directory `meta`. 
There, you find a list of links as well as a download script for automatition.
As the meta-data is provided as an xml file per sample, we also implemented a tool `xml_to_csv.py`
to convert all these xml files into a single csv table, which is provided in `data.csv`.
Lastly, `data_simple.csv` contains a selection of the meta-data columns which we used for the paper.

Additional Files
-------------------------

`histogram.ods`: length histroams of the sequences. Used for deciding which sequence lengths to exclude.
