#!/bin/bash

genesis/bin/apps/label_matrix data/hmp/meta/label_body_site.txt data/hmp/meta/label_body_site_list.txt 01_backbone/10_kmeans/hmp_Bact_Constr_backbone/emd_assignments.csv 01_backbone/10_kmeans/hmp_Bact_Constr_backbone/emd_body_site_orrd.bmp


genesis/bin/apps/label_matrix data/hmp/meta/label_body_site.txt data/hmp/meta/label_body_site_list.txt 01_backbone/10_kmeans/hmp_Bact_Constr_backbone/imbalance_assignments.csv 01_backbone/10_kmeans/hmp_Bact_Constr_backbone/imb_body_site_orrd.bmp

pixel2svg-0.3.0/pixel2svg.py --overlap 01_backbone/10_kmeans/hmp_Bact_Constr_backbone/emd_constr_body_site_orrd.bmp
pixel2svg-0.3.0/pixel2svg.py --overlap 01_backbone/10_kmeans/hmp_Bact_Constr_backbone/imb_consr_body_site_orrd.bmp
