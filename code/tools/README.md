Overview
-------------------------

These programs are tools that might be reusable for other purposes as well.

`fasta_cleanup`
-------------------------

Remove all non ACGT chars from the sequences of a fasta file, and replace U by T.

`fasta_labels`
-------------------------

Write out all sequence names of a fasta file.

`fasta_to_phylip`
-------------------------

Simple converter. Sequences have to be aligned (i.e., all have the same length),
because phylip...

`fasta_to_uppper`
-------------------------

Make sequence all upper case.

`fasta_unalign`
-------------------------

Remove all gaps from a fasta file, that is, "unalign" it.

`jplace_to_bplace`
-------------------------

Converter into our internal `bplace` format, which is needed for some of the programs.

`jplace_tree_info`
-------------------------

Print information about the reference tree of a jplace file,
like average branch length etc.

`mat_to_bmp`
-------------------------

Visualize a matrix (e.g., a pairwise distance matrix) as a grey color heatmap.

`msa-[head|tail]`
-------------------------

Get the first or last `n` sequences of a fasta file.

`merge_jplace_files`
-------------------------

Take a list of jplace files and merge them (required them to have the same ref tree).
Write to the file that is the last command line argument.

`newick_compare`
-------------------------

Compare two newick files and visualize their confilcting branches.
This is baiscally a viz of the RF distances.

`newick_info`
-------------------------

Print information about the reference tree of a newick file,
like average branch length etc.

`phylip_to_fasta`
-------------------------

Simple converter.
