#!/bin/bash

# Perform a full simulation for a given generation method, number of
# replicates, number of injected simulated relationships, and 
# simulation plan.

GENERATION_METHOD=$1
REPLICATES=$2
N_RELATIONS=$3
SIM_PLAN=$4

export REPLICATES
export GENERATION_METHOD

pushd ${PROJ_DIR}

PLAN_TAG=$(basename $SIM_PLAN .csv)

./sh/conserved-roc.sh $N_RELATIONS $SIM_PLAN > work/${GENERATION_METHOD}-${REPLICATES}-${N_RELATIONS}-${PLAN_TAG}.out 2>&1

popd

