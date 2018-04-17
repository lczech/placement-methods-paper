#!/bin/bash

BASEDIR=hmp_kmeans_test
TESTDIR=hmp_kmeans_test/kmeans_elbow

/usr/bin/time -v ${TESTDIR}/genesis/bin/apps/kmeans_jplace_elbow 40 32 ${BASEDIR}/samples_jplace ${TESTDIR} &> log.txt

