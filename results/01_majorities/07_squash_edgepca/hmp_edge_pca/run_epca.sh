#!/bin/bash

BASEDIR=hmp_edge_pca

/usr/bin/time -v ${BASEDIR}/genesis/bin/apps/bplace_epca ${BASEDIR}/samples_bplace ${BASEDIR} &> log_epca.txt

