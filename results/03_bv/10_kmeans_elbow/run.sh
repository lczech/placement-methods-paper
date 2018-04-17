#!/bin/bash

BASEDIR=bv_kmeans_test
TESTDIR=bv_kmeans_test/kmeans_elbow

/usr/bin/time -v ${TESTDIR}/genesis/bin/apps/kmeans_jplace_elbow 10 16 ${BASEDIR}/orig_queries_jplace_clean ${TESTDIR} &> log.txt

