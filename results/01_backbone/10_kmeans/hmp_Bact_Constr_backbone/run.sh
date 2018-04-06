#!/bin/bash

genesis/bin/apps/kmeans_bplace 12 19 06_samples/hmp_Bact_Constr_backbone/samples_bplace 10_kmeans/hmp_Bact_Constr_backbone &> log.txt

mkdir figs
mv *.nexus *.svg figs/

