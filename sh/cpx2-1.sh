#!/bin/bash
# Run CPX2 on one set of trajectories.

ANA_DIR=${PROJ_DIR}/work/R1

T1=${ANA_DIR}/FBgn0030774-FBgn0029937-R1.trj
T2=${ANA_DIR}/FBgn0030774-FBgn0029937-__X__R1.trj

echo Trajectory 1:
cat $T1

echo Trajectory 2:
cat $T2

${PROJ_DIR}/bin/CPX2-linux -M comparison -p 1 -g 1 -K 0 -J 0  -1 $T1 -2 $T2

