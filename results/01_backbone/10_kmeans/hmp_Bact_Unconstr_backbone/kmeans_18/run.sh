#!/bin/bash

BASEDIR=hmp_kmeans_test

/usr/bin/time -v ${BASEDIR}/genesis/bin/apps/kmeans_bplace 16 18 ${BASEDIR}/samples_bplace ${BASEDIR}/kmeans_18 &> log.txt

