# Utility functions.

util.util <- function() {
  # I exist only as a flag.
}

# For getting area under rough ROC curves.
require(pracma)

PROJ_DIR = Sys.getenv('PROJ_DIR')
if (PROJ_DIR == "") {
    cat("Environment variable PROJ_DIR should be set to the project root directory.\n")
    quit('no')
}

aucFromTable <- function(roc.dat) {
    max.fp.proportion = roc.dat[nrow(roc.dat),'false.proportion']
    AUC = trapz(roc.dat[,'false.proportion'],roc.dat[,'true.proportion'])/max.fp.proportion
    return(AUC)
}

plotRocFromTable <- function(roc.dat) {
    AUC = aucFromTable(roc.dat)
    plot(x=roc.dat[,'false.proportion'],y=roc.dat[,'true.proportion'],xlab="False positive proportion",ylab="True positive proportion")
    title(main=paste("AUROC ~ ",AUC,sep='',collapse=''))
    lines(x=roc.dat[,'false.proportion'],y=roc.dat[,'true.proportion'],type="l")
}

# Convert a Flybase ID to a gene name, if possible. Otherwise
# return the original Flybase ID.
gnameFile <- 'data/Dutta_FlyBase_IDs.txt'
gnameTab <- read.csv(gnameFile,sep='\t',stringsAsFactors=FALSE)

gname <- function(fbid) {
    nm = gnameTab[gnameTab[,2] == fbid,4]
    if (any(is.na(nm)) || any(is.null(nm)) || any(length(nm)==0)) {
        nm = fbid
    }
    return(nm)
}

fdr_BH <- function(data,BH.alpha=0.25) {
    sorted.data = data[order(data[,'p.value']),]
    sorted.data[['BH.significant']] <- rep(FALSE,nrow(sorted.data))
    for (jj in 1:nrow(sorted.data)) {
        if (sorted.data[jj,'p.value'] < (jj*BH.alpha/nrow(sorted.data))) {
            sorted.data[jj,'BH.significant'] <- TRUE
        }
    }
    return(sorted.data)
}

pathAndName <- function(path) {
    parts = strsplit(path,'/',fixed=T)[[1]]
    last = parts[length(parts)]
    rest = paste(parts[-length(parts)],collapse='/')
    return(c(rest,last))
}

emptyRow = function(nReplicates=2) {
    nr = 0
    ll = data.frame()
    ll[["FBgn"]] = character(0)
    ll[["CG number"]] = character(0)
    ll[["Gene ID"]] = character(0)
    for (cellType in c("EC","EB","EE","ISC")) {
        for (region in c("R1","R2","R3","R4","R5")) {
            for (rep in 1:nReplicates) {
                cname = paste(cellType,"_",region,"_",as.character(rep),collapse="",sep="")
                ll[[cname]] = numeric(0)
            }
        }
    }
    return(ll)
}

