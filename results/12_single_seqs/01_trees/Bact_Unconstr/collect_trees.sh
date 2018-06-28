#!/bin/bash

LH_STRING="Final GAMMA-based Score of best tree "

LH_FILE="LH_bests"
SORTED_LH_FILE="LH_bests_sorted"
rm -f ${LH_FILE} ${SORTED_LH_FILE}

for dir in `ls -d s*`
do
    if [ -f "$dir" ]
    then
      continue
    fi

    echo $dir
    cd $dir

    INFO_FILE=`ls RAxML_info.*`
    TREE_FILE=`ls RAxML_bestTree.*`

    if [ -z "$INFO_FILE" ]
    then
      echo "Skipping $dir"
      cd ..
      continue
    fi

    LH=`grep "${LH_STRING}" ${INFO_FILE} | sed s/"${LH_STRING}"//g`
    echo "$LH ${dir}/${TREE_FILE}" >> ../$LH_FILE

    cd ..
done

sort -n -r $LH_FILE > ${SORTED_LH_FILE}
BEST_RESULT=`head -n 1 ${SORTED_LH_FILE} | tr -s " " "\n" | sed -n '2p'`
echo "Best tree: $BEST_RESULT"
cp ${BEST_RESULT} best_tree.newick
