#!/bin/bash

BASEDIR="path/to/here"

echo "Start at `date`"
echo

/usr/bin/time -v ${BASEDIR}/pplacer/guppy kr --out-dir ${BASEDIR} -o "guppy_kr.txt" ${BASEDIR}/orig_queries_jplace_clean/*

echo
echo "End at `date`"
