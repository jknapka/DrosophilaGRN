#!/bin/bash

IN_FILE=$1
OUT_FILE=$2

# Find interactions in the physical_interactions file whose
# source and target genes both appear in the Dutta dataset.
lines=$(egrep "^FBgn" ${IN_FILE})
while read -r line ; do
    g1=$(echo $line | awk '{print $1}');
    g2=$(echo $line | awk '{print $3}');
    if grep -q $g1 data/Table-S1.csv && grep -q $g2 data/Table-S1.csv ; then
        echo $g1,$g2 ;
    fi
done <<< "$lines" > ${OUT_FILE}

