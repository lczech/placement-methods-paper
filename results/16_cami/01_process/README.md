Overview
================

### preprocessing

call `do_all.sh` from the `19122017_mousegut_scaffolds` folder, specifiying all sample paths:
```
../place/scripts/do_all.sh 2017.12.29_11.37.26_sample_*
```
This will preprocess each sample into a <sampledir>/reads/sample_x.16S_rRNA and <same>...18S_rRNA that includes only the ssu reads

Next, we use `gappa` to globally de-replicate the data (now we are operating from the `place` folder):
```
mkdir -p query/abundance/
gappa prepare chunkify --abundances-out-dir query/abundance/ --chunks-out-dir query/ --chunk-size 1000000 --fasta-path 2017.12.29_11.37.26_sample_*/reads/sample_*.16S_rRNA
```

### alignment

for unconstrained ref tree:
```
papara -t tree/bact_ref_pkg/majorities_Bacteria_Unconstr_best_tree.newick -s tree/bact_ref_pkg/Bacteria_sequences.fasta.reduced -q query/chunk_0.fasta -n uncons_aligned -r
```
and constrained:
```
papara -t tree/bact_ref_pkg/majorities_Bacteria_Constr_best_tree.newick -s tree/bact_ref_pkg/Bacteria_sequences.fasta.reduced -q query/chunk_0.fasta -n cons_aligned -r
```
then prepare for placement by converting and splitting off the queries. Again first unconstrained:
```
cd query
epa-ng --split ../tree/bact_ref_pkg/Bacteria_sequences.fasta.reduced papara_alignment.uncons_aligned
mv query.fasta uncons_query.fasta
```
then constrained:
```
epa-ng --split ../tree/bact_ref_pkg/Bacteria_sequences.fasta.reduced papara_alignment.cons_aligned
mv query.fasta cons_query.fasta
cd -
```

### placement

unconstrained:
```
mkdir -p jplace/uncons/
epa-ng --tree tree/bact_ref_pkg/majorities_Bacteria_Unconstr_best_tree.newick --ref-msa tree/bact_ref_pkg/Bacteria_sequences.fasta.reduced.fasta --query query/uncons_query.fasta --model tree/bact_ref_pkg/info/RAxML_info.info --outdir jplace/uncons
```
constrained:
```
mkdir -p jplace/uncons/
epa-ng --tree tree/bact_ref_pkg/majorities_Bacteria_Constr_best_tree.newick --ref-msa tree/bact_ref_pkg/Bacteria_sequences.fasta.reduced.fasta --query query/cons_query.fasta --model tree/bact_ref_pkg/info/RAxML_info.info --outdir jplace/cons
```

de-dereplication to obtain per-sample result jplace files
```
cd jplace/uncons # OR cd jplace/cons depending which you want
mkdir -p samples
gappa prepare unchunkify --abundances-path ../../query/abundance/ --jplace-path epa_result.jplace --out-dir samples
cd -
```

### taxonomic assignment

(from the `place` folder)
unconstrained:
```
./assign_all.sh uncons
```
constrained:
```
./assign_all.sh cons
```
