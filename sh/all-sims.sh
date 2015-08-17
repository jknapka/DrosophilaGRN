#!/bin/bash
#
# Run all simulations of interest: gen, rgen, cval for 2, 3, 4, 5, 10 replicates,
# injecting simDiffEB-clear.csv, simDiffEB-noisy-03, and simDiffEB-noisy-05.csv.
pushd ${PROJ_DIR}

for gmeth in gen rgen cval ; do
    for repls in 2 3 4 5 10 ; do

        for simPlan in simDiffEB-clear.csv simDiffEB-noisy-03.csv simDiffEB-noisy-05.csv ; do

            echo Running $gmeth $simPlan with $repls replicates.

            ./sh/do-simulation.sh ${gmeth} ${repls} 10 ${simPlan}

        done

    done
done

popd
