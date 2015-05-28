#!/bin/bash
# Run CPX2 on one set of trajectories.

ANA_DIR=${PROJ_DIR}/work/$2

T1=${ANA_DIR}/$1-$2.trj
T2=${ANA_DIR}/$1-__X__$2.trj

echo Trajectory 1:
cat $T1

echo Trajectory 2:
cat $T2

${PROJ_DIR}/bin/glnsp -M comparison -p 1 -g 1 -K 0 -J 0  -1 $T1 -2 $T2

