#!/bin/bash

# Build a randomized data set containing some
# simulated interactions and evaluate the effects
# of varying the Benjamini-Hochberg FDR parameter.

if [ "" == "${PROJ_DIR}" ] ; then
    echo PROJ_DIR is not set.
    exit 1
fi

cd $PROJ_DIR

# Build the random data set.
if [ "" == "${REPLICATES}" ] ; then REPLICATES=2 ; fi
if [ "" == "${GENERATION_METHOD}" ] ; then GENERATION_METHOD=gen ; fi 


JOB_DIR=${PROJ_DIR}/work/conserved-ROC/${GENERATION_METHOD}-${REPLICATES}-$(date +%s)
mkdir -p $JOB_DIR

echo GENERATION METHOD: $GENERATION_METHOD > ${JOB_DIR}/params.txt
echo REPLICATES: $REPLICATES >> ${JOB_DIR}/params.txt
echo PLANS: $@ >> ${JOB_DIR}/params.txt

Rscript --no-restore ${PROJ_DIR}/r/randomGeneRow.r $GENERATION_METHOD ${JOB_DIR}/random.csv $REPLICATES
#Rscript --no-restore ${PROJ_DIR}/r/randomGeneRow.r cval ${JOB_DIR}/random.csv $REPLICATES

if [ ! "" == "$1" ] ; then
    CMD_ARGS=()
    while [ ! "" == "$1" ] ; do
        N_INTERACTIONS=$1
        INTERACTION_FILE=$2
        echo Injecting $N_INTERACTIONS instances of $INTERACTION_FILE
        # Inject some simulated interactions.
        for (( x=0 ; x<= $N_INTERACTIONS ; ++x )) ; do
            CMD_ARGS+=( ${INTERACTION_FILE} )
        done
        shift 2
    done
    echo Running the injection.
    Rscript --no-restore ${PROJ_DIR}/r/randomGeneRow.r inject ${JOB_DIR}/random.csv ${JOB_DIR}/injected.csv "${CMD_ARGS[@]}"
else
    # Just check that we get the expected results with NO simulated interactions.
    cp ${JOB_DIR}/random.csv ${JOB_DIR}/injected.csv
    touch ${JOB_DIR}/simint-injected.csv
fi

# Evaluate conserved interactions. This also discretizes the
# dataset and writes a file we use to do differential analysis.
Rscript --no-restore ${PROJ_DIR}/r/conserved-interactions.r ${JOB_DIR}/injected.csv ${JOB_DIR}/simint-injected.csv > ${JOB_DIR}/conserved.out
mv Rplots.pdf ${JOB_DIR}/conserved-roc.pdf

${PROJ_DIR}/sh/differential-roc.sh ${JOB_DIR}

echo Done.

