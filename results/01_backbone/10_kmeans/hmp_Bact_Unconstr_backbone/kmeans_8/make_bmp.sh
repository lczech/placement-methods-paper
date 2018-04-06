#!/bin/bash

#genesis/bin/apps/label_matrix data/hmp/meta/label_body_site.txt data/hmp/meta/label_body_site_list_8.txt 01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8/emd_assignments.csv 01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8/emd_body_site_orrd.bmp

#genesis/bin/apps/label_matrix data/hmp/meta/label_body_site.txt data/hmp/meta/label_body_site_list_8.txt 01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8/imbalance_assignments.csv 01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8/imb_body_site_orrd.bmp

#pixel2svg-0.3.0/pixel2svg.py --overlap 01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8/emd_body_site_orrd.bmp
#pixel2svg-0.3.0/pixel2svg.py --overlap 01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8/imb_body_site_orrd.bmp

rm -r *_8.svg *_8.bmp

data_dir=01_backbone/10_kmeans/hmp_Bact_Unconstr_backbone/kmeans_8

label_body_site=data/hmp/meta/label_body_site.txt
label_body_site_list=data/hmp/meta/label_body_site_list_8.txt
emd_assignments=${data_dir}/emd_assignments.csv
imbalance_assignments=${data_dir}/imbalance_assignments.csv

genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${emd_assignments} ${data_dir}/emd_body_site_orrd_counts_8.bmp
genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${imbalance_assignments} ${data_dir}/imb_body_site_orrd_counts_8.bmp

genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${emd_assignments} ${data_dir}/emd_body_site_orrd_norm_8.bmp normalize
genesis/bin/apps/label_matrix ${label_body_site} ${label_body_site_list} ${imbalance_assignments} ${data_dir}/imb_body_site_orrd_norm_8.bmp normalize

pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/emd_body_site_orrd_counts_8.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/imb_body_site_orrd_counts_8.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/emd_body_site_orrd_norm_8.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap ${data_dir}/imb_body_site_orrd_norm_8.bmp
