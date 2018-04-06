#!/bin/bash

for line in `cat ../00_data/sample_names`; do
	grep "${line}	" meta_raw >> meta_all
done
