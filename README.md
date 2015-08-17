# DrosophilaGRN
Code for analyzing Drosophila transcriptome data

You will need to build a recent version of GLN (executable name glnsp) and
place glnsp in ./bin/

You will need to create ./work/ before doing anything significant.

Directories:

r/  contains R code to do parts of the simulation and analysis.

sh/  contains driver scripts that invoke glnsp to perform differntial interaction analysis,
and also scripts to perform various simulations. sh/all-sims.sh runs a large number of
simulations with different parameters. Results end up in work/conserved-ROC (even though
much of the analysis is actuall concerned with detecting differential interactions).

simulation/  contains obsolete Python code to build simulated data sets. It is
all superseded by r/randomGeneRow.r


