if (!exists("util.util",mode="function")) {
    # print("Reading grn_util.r")
    source("r/grn_util.r")
}

args <- commandArgs(trailingOnly=TRUE)

if (args[1] == '--from-table') {
    cat("rocPlot.r reading data file",args[2],'\n')
    flush.console()
    dat = read.csv(args[2],sep=',',header=T,stringsAsFactors=F)
    cat("rocPlot.r reading interaction file",args[3],'\n')
    flush.console()
    simints = read.csv(args[3],sep=',',header=F,stringsAsFactors=F,col.names=c('FBid1','FBid2'))
    simdat = addSimulationFlag(dat,simints,fbid.col.1=2,fbid.col.2=3)
    roc.list = buildROC(simdat)
    roc.dat = roc.list[['ROC']]

    print(roc.dat)

    # Write out the interactions where BH.alpha == 0.2
    interactions = roc.list[['interactions']]
    point2flags = unlist(lapply(interactions,FUN=function(x){return(x[1,'BH.alpha']>= 0.2)}))
    point2idxs = which(point2flags == T)
    cat("point2flags",point2flags[1:10],'\npoint2idxs',point2idxs[1:10],'\n')
    fdr2.interactions = interactions[[min(point2idxs)]]
    cat("Printing interactions where BH.alpha =",fdr2.interactions[1,'BH.alpha'],'\n')
    flush.console()
    write.table(fdr2.interactions,args[2],sep=',',col.names=T,row.names=F,quote=F)
} else {
    inFile = args[1]
    roc.dat = read.csv(inFile,sep=',',header=T,stringsAsFactors=F)
}

plotRocFromTable(roc.dat)

