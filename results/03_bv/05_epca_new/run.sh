#!/bin/bash

echo "Start at `date`"

/usr/bin/time -v genesis/bin/apps/jplace_epca_vis 03_bv/03_epa/orig_queries_jplace_clean 03_bv/05_epca_new

echo "End at `date`"
