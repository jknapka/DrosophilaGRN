#!/bin/bash

DATA_TEMPLATE_FILE=$1
shift
NSIMS=$1
shift
DEST_FILE=$1
shift
REPLICATES=$1
shift

SIM_ARGS="$@ -r ${REPLICATES}"

WORK_DIR=$(mktemp -d ./sim-work-dir-XXXXXXXXXX)

CUR_FILE=${DATA_TEMPLATE_FILE}
declare -a INJECT_FILES

if [ ! -e ${DATA_TEMPLATE_FILE} ] ; then
    Rscript --no-restore ${PROJ_DIR}/r/randomGeneRow.r gen ${DATA_TEMPLATE_FILE} ${REPLICATES}
fi

for (( x=1 ; x<=$NSIMS ; ++x )) ; do

    IFILE=${WORK_DIR}/i${x}.csv
    python ${PROJ_DIR}/simulation/makeDifferentialRow.py $SIM_ARGS > ${IFILE}
    INJECT_FILES[${x}]=${IFILE}

done

egrep -o "^FBgn[0-9]+" ${WORK_DIR}/i* | cut -d: -f2 | xargs -n 2 echo | tr \  , > ${WORK_DIR}/interactions.txt
Rscript --no-restore ${PROJ_DIR}/r/randomGeneRow.r inject ${DATA_TEMPLATE_FILE} ${DEST_FILE} ${INJECT_FILES[@]}

