#!/bin/bash

PROGRAM=genesis/bin/apps/dispersion_trees
BASE_DIR="/path/to/08_emd"

for DATASET in bv hmp tara; do
  for TREE in Bact Euks Gen; do
    for CONSTR in Constr Unconstr; do

      WORK_DIR=${BASE_DIR}/${DATASET}_${TREE}_${CONSTR}_backbone
      if [ ! -d ${WORK_DIR} ]; then
        continue
      fi

      TREE_FILE=${WORK_DIR}/avg_tree.newick
      EDGE_FILE=${WORK_DIR}/edge_weights.mat
      IMB_FILE=${WORK_DIR}/imbalance.mat
      OUT_PATH=${WORK_DIR}/

      echo ${WORK_DIR}
      echo "==========================="
      ${PROGRAM} ${TREE_FILE} ${EDGE_FILE} ${IMB_FILE} ${OUT_PATH}

    done
  done
done

exit

# find . -name disp_edge_* | xargs rm -v

#./dispersion_trees 08_emd/bv_Gen_Constr_backbone/avg_tree.newick 08_emd/bv_Gen_Constr_backbone/edge_weights.mat 08_emd/bv_Gen_Constr_backbone/imbalance.mat 08_emd/bv_Gen_Constr_backbone/
