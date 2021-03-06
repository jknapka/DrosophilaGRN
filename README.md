# DrosophilaGRN
Code for analyzing Drosophila transcriptome data and simulating expression data for
regulatory interactions between genes.

All code in this tree expects the environment variable PROJ_DIR to
be set to the full path of the directory containing this README.md file.

There is quite a bit of code in here that assumes that expression data
is available in the dataset for cell types named "EB", "EC", "EE", and
"ISC". In retrospect, this was unwise! The solution probably involves
selecting data columns of interest using a regular expression that can be
customized for any dataset.

There is no real data for analysis in this repository, because the data
I analyzed had not been published at the time I was doing this work. To
be useful, you will need to create a ./data directory and put some things
there.

Some of this code assumes that the "real" data lives in a table called
./data/Table-S1.csv. Most things can take an argument to specify the table
to process, however.

You need a file called ./data/interactions-present-in-data.txt, which is
really a comma-separated file with two fields per line. Each line names
two genes that interact.

The code for converting gene IDs to names in r/grn_util.r depends on some
data that is not present here. It probably isn't really necessary and
references to the gname() function should be removed.

You will need to build a recent version of GLN (executable name glnsp) and
place glnsp in ./bin/

You will need to create ./work/ before doing anything significant.

Directories:

r/  contains R code to do parts of the simulation and analysis.

sh/  contains driver scripts that invoke glnsp to perform differntial interaction analysis,
and also scripts to perform various simulations. sh/all-sims.sh runs a large number of
simulations with different parameters. Results end up in work/conserved-ROC (even though
much of the analysis is actually concerned with detecting differential interactions).

simulation/  contains obsolete Python code to build simulated data sets. That code is
all superseded by r/randomGeneRow.r

simulation/plans contains simulation plans that can be used with
"r/randomGeneRow.r inject" to add known interactions to a randomized dataset.
The comment for the injectData2 function in r/randomGeneRow.r explains the format 
of a simulation plan file.  Briefly,here is a simulation plan which will be
explained below:

```
celltype,    level,    maxlevel,    tlevel,    maxtlevel,    noise
EB,          1,        3,           1,         10,           0.2
EB,          2,        3,           2,         10,           0.2
EB,          3,        3,           10,        10,           0.3
OTHER,       1,        1,           1,         2,            1.0
OTHER,       2,        2,           1,         2,            1.0
OTHER,       3,        3,           1,         2,            1.0
OTHER,       4,        4,           1,         2,            1.0
OTHER,       5,        5,           1,         2,            1.0
```

The columns are the cell type, source gene level, maximum source gene level,
target level to use when the interaction is simulated for the corresponding
source level, maximum target level, and "house noise" level.

To begin implementing a simulation plan, source and target gene expression
profiles are chosen at random from the genes present in the original
"real" dataset. (The expression profiles are written out by
conserved-interactions.r to a file called data/geneStats.csv. That
has to happen before a simulation can be run.) Then source and target
gene IDs are chosen from the simulated dataset; the data for these genes
will be overwritten by the simulated expression values.

When simulating expression data for a cell of type EB, a random row from
the celltype=EB rows will be chosen. The source gene's level will be set
to the level value (mapped to the actual range of expression levels
of the selected target gene by interpolation from the range [1..maxlevel]).
The target gene's level will be set to the corresponding tlevel value
(interpolated from the range [1..maxtlevel], and the "house noise" model
will be applied to possibly change the target level according to the
value in the noise column. The functional relationship indicated by this
plan for EB cells is:
```
1/3 -> 1/10          0.333 -> 0.1
2/3 -> 3/10   (or)   0.667 -> 0.3
3/3 -> 10/10         1.0   -> 1.0
```
which will subsequently be linearly interpolated to the range of levels 
actually discovered by Ckmeans for the source and target genes in the
"real" dataset.

When simulating expression data for any other cell type, the OTHER row
is used. A source level is chosen at random, and the target level is
always 1. However, since the noise value is 1.0, the effect is that a
level will be chosen uniformly at random from the target gene's range
of quantized levels. Thus, this plan simulates an interaction with a
strong exponential-shaped signal for EB cells, and random expression data
for all other cell types.
