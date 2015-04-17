#!/bin/bash
lines=$(egrep "^FBgn" data/physical_interactions_fb_2015_01.tsv)
while read -r line ; do
    g1=$(echo $line | awk '{print $1}');
    g2=$(echo $line | awk '{print $3}');
    if grep -q $g1 data/Table-S1.csv && grep -q $g2 data/Table-S1.csv ; then
        echo $g1,$g2 ;
    fi
done <<< "$lines" > interactions-present-in-data.txt

