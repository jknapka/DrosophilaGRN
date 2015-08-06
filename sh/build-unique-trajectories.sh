#!/bin/bash

INFILE=$1

cd $PROJ_DIR

for category in R1 R2 R3 R4 R5 EC EE EB ISC ; do
    Rscript --no-restore r/unique-interactions.r $category $INFILE
done
