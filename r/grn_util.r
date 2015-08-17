# Utility functions.

util.util <- function() {
  # I exist only as a flag.
}

# For getting area under rough ROC curves.
require(pracma)

# Ensure PROJ_DIR is set appropriately.
PROJ_DIR = Sys.getenv('PROJ_DIR')
if (PROJ_DIR == "") {
    cat("Environment variable PROJ_DIR should be set to the project root directory.\n")
    quit('no')
}

# Compute the area under a ROC curve from the corresponding values
# of false.proportion and true.proportion in the given data.frame.
aucFromTable <- function(roc.dat) {
    max.fp.proportion = roc.dat[nrow(roc.dat),'false.proportion']
    if (max.fp.proportion == 0) { max.fp.proportion = 1.0 }
    AUC = trapz(roc.dat[,'false.proportion'],roc.dat[,'true.proportion'])/max.fp.proportion
    return(AUC)
}

# Plot the ROC from a data.frame having true.proportion and false.proportion
# columns, indicating the proportion of true and false positives.
plotRocFromTable <- function(roc.dat) {
    AUC = aucFromTable(roc.dat)
    plot(x=roc.dat[,'false.proportion'],y=roc.dat[,'true.proportion'],xlab="False positive proportion",ylab="True positive proportion")
    title(main=paste("AUROC ~ ",AUC,sep='',collapse=''))
    lines(x=roc.dat[,'false.proportion'],y=roc.dat[,'true.proportion'],type="l")
}

# Convert a Flybase ID to a gene name, if possible. Otherwise
# return the original Flybase ID.
gnameFile <- file.path(PROJ_DIR,'data/Dutta_FlyBase_IDs.txt')
gnameTab <- read.csv(gnameFile,sep='\t',stringsAsFactors=FALSE)

gname <- function(fbid) {
    nm = gnameTab[gnameTab[,2] == fbid,4]
    if (any(is.na(nm)) || any(is.null(nm)) || any(length(nm)==0)) {
        nm = fbid
    }
    return(nm)
}

# Add Benjamini-Hochberg FDR indication as a new column in the given
# data.frame, which must have a p.value column.
fdr_BH <- function(data,BH.alpha=0.25) {
    #cat("Performing FDR with BH.alpha=",BH.alpha,"and p.values=",data[['p.value']])
    sorted.data = data[order(data[,'p.value']),]
    sorted.data[['BH.significant']] <- rep(FALSE,nrow(sorted.data))
    for (jj in 1:nrow(sorted.data)) {
        if ((BH.alpha == 1.0) || (sorted.data[jj,'p.value'] < (jj*BH.alpha/nrow(sorted.data)))) {
            sorted.data[jj,'BH.significant'] <- TRUE
        }
    }
    return(sorted.data)
}

# Compute the path and base name of a file.
pathAndName <- function(path) {
    parts = strsplit(path,'/',fixed=T)[[1]]
    last = parts[length(parts)]
    rest = paste(parts[-length(parts)],collapse='/')
    return(c(rest,last))
}

# Build an empty data row for the given number of replicates.
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

# Build the ROC table for a given set of observations.
buildROC <- function(interactions,BH.resolution=0.01) {
    roc.data = data.frame(FDR=numeric(0),true.positives=numeric(0),true.proportion=numeric(0),false.positives=numeric(0),false.proportion=numeric(0))
    allInteractionSets = list()
    for (BH.alpha in seq(0.0,1.0,BH.resolution)) {
        sorted.interactions = fdr_BH(interactions,BH.alpha=BH.alpha)

        sorted.interactions[['n']] = 1:nrow(sorted.interactions)
        sorted.interactions[['BH.alpha']] = rep(BH.alpha,nrow(sorted.interactions))

        allInteractionSets[[length(allInteractionSets)+1]] = sorted.interactions

        #print("Adjusted interaction significance:")
        #print(sorted.interactions)

        true.positives = nrow(sorted.interactions[sorted.interactions[,'BH.significant'] & sorted.interactions[,'simulated'],])
        n.true = nrow(sorted.interactions[sorted.interactions[,'simulated'],])
        false.positives = nrow(sorted.interactions[sorted.interactions[,'BH.significant'] & !sorted.interactions[,'simulated'],])
        n.false = nrow(sorted.interactions)-n.true

        cat("  At FDR",BH.alpha,"found",true.positives,"/",n.true,"true positives and",false.positives,"/",n.false,"false positives\n")

        pTrue = 0.0
        pFalse = 0.0
        if (n.true > 0) {
            pTrue = true.positives/n.true
        }
        if (n.false > 0) {
            pFalse = false.positives/n.false
        }

        #cat("lengths: BH.alpha",length(BH.alpha),"true.positives",length(true.positives),"true.proportion",length(pTrue),
        #    "false.positives",length(false.positives),"false.proportion",length(pFalse))
        roc.row = list(FDR=BH.alpha,true.positives=true.positives,true.proportion=pTrue,false.positives=false.positives,false.proportion=pFalse)

        #print("trying to rbind")
        #print(roc.data)
        #print("   to")
        #print(roc.row)
        #flush.console()

        roc.data = rbind(roc.data,roc.row)
    }
    return(list(ROC=roc.data,interactions=allInteractionSets))
}

# Add a column indicating which interactions are simulated.
addSimulationFlag <- function(interactions,simulated_interactions,fbid.col.1=1,fbid.col.2=2) {
    cat("Adding simulation flag to",nrow(interactions),"interactions.\n")
    interactions[['simulated']] <- rep(FALSE,nrow(interactions))

    if (nrow(simulated_interactions) < 1) return(interactions)
    if (is.null(simulated_interactions)) return(interactions)

    for (jj in 1:nrow(interactions)) {
        for (kk in 1:nrow(simulated_interactions)) {
            # Genes in this interation.
            ig1 = interactions[jj,fbid.col.1]
            ig2 = interactions[jj,fbid.col.2]
            # Genes in this simulated interaction.
            sg1 = simulated_interactions[kk,1]
            sg2 = simulated_interactions[kk,2]
            if ((ig1==sg1) && (ig2==sg2)) {
                interactions[jj,'simulated'] <- TRUE
            }
        }
    }
    return(interactions)
}

