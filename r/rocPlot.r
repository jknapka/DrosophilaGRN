if (!exists("util.util",mode="function")) {
    # print("Reading grn_util.r")
    source("r/grn_util.r")
}

args <- commandArgs(trailingOnly=TRUE)

inFile = args[1]

roc.dat = read.csv(inFile,sep=',',header=T)
plotRocFromTable(roc.dat)

