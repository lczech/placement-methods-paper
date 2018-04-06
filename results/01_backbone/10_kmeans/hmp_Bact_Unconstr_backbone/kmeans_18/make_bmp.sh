#!/bin/bash

rm -r *.svg *.bmp

data_dir=01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_18

label_body_site=data/hmp/meta/label_body_site.txt
label_body_site_list=data/hmp/meta/label_body_site_list_18.txt
emd_assignments=${data_dir}/emd_assignments.csv
imbalance_assignments=${data_dir}/imbalance_assignments.csv

genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${emd_assignments} ${data_dir}/emd_body_site_orrd_counts_18.bmp
genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${imbalance_assignments} ${data_dir}/imb_body_site_orrd_counts_18.bmp

genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${emd_assignments} ${data_dir}/emd_body_site_orrd_norm_18.bmp normalize
genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${imbalance_assignments} ${data_dir}/imb_body_site_orrd_norm_18.bmp normalize

pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/emd_body_site_orrd_counts_18.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/imb_body_site_orrd_counts_18.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/emd_body_site_orrd_norm_18.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/imb_body_site_orrd_norm_18.bmp
