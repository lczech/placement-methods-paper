#!/bin/bash

genesis/bin/apps/kmeans_bplace 6 3 06_samples/bv_Bact_Constr_backbone/samples_bplace 10_kmeans/bv_Bact_Constr_backbone &> log.txt

mkdir figs
mv *.nexus *.svg figs/

