#!/bin/bash

LABEL_MATRIX=genesis/bin/apps/label_matrix
KMEANS_DIR="10_kmeans"
META_DIR="data"

#####################################################
#    HMP
#####################################################

DATASET="hmp"
for TREE in Bact Gen; do
  for CONSTR in Constr Unconstr; do
    for META in "body_site"; do
      for DIST in emd imbalance; do
        echo ${TREE} ${CONSTR} ${META}

        LABEL_FILE=${META_DIR}/${DATASET}/meta/label_${META}.txt
        LABEL_LIST=${META_DIR}/${DATASET}/meta/label_${META}_list.txt
        ASSIGN_FILE=${KMEANS_DIR}/${DATASET}_${TREE}_${CONSTR}_backbone/${DIST}_assignments.csv
        OUT_FILE=${KMEANS_DIR}/${DATASET}_${TREE}_${CONSTR}_backbone/${DIST}_${META}_viridis.bmp

        ${LABEL_MATRIX} ${LABEL_FILE} ${LABEL_LIST} ${ASSIGN_FILE} ${OUT_FILE}
      done
    done
  done
done

exit

#####################################################
#    BV
#####################################################

DATASET="bv"
for TREE in Bact Gen; do
  for CONSTR in Constr Unconstr; do
    for META in amsel nugent; do
      for DIST in emd imbalance; do
        echo ${TREE} ${CONSTR} ${META}

        LABEL_FILE=${META_DIR}/${DATASET}/meta/label_${META}.txt
        LABEL_LIST=${META_DIR}/${DATASET}/meta/label_${META}_list.txt
        ASSIGN_FILE=${KMEANS_DIR}/${DATASET}_${TREE}_${CONSTR}_backbone/${DIST}_assignments.csv
        OUT_FILE=${KMEANS_DIR}/${DATASET}_${TREE}_${CONSTR}_backbone/${DIST}_${META}_viridis.bmp

        ${LABEL_MATRIX} ${LABEL_FILE} ${LABEL_LIST} ${ASSIGN_FILE} ${OUT_FILE}
      done
    done
  done
done

#./label_matrix data/bacterial_vaginosis/label_nugent.txt 10_kmeans/bv_Bact_Unconstr_backbone/imbalance_assignments.csv bv_Bact_Unconstr_backbone.bmp

