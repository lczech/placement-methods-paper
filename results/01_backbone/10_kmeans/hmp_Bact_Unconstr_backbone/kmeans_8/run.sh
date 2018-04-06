#!/bin/bash

BASEDIR=hmp_kmeans_test

/usr/bin/time -v ${BASEDIR}/genesis/bin/apps/kmeans_bplace 16 8 ${BASEDIR}/samples_bplace ${BASEDIR}/kmeans_8 &> log.txt

