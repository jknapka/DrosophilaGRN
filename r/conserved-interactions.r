require('Ckmeans.1d.dp')
require('FunChisq')

if (!exists("util.util",mode="function")) {
    print("Reading grn_util.r")
    source("r/grn_util.r")
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
    inFile <- args[1]
    gsFile <- paste(inFile,"-geneStats.csv",sep="",collapse="")
    if (length(args) > 1) {
        simulated_interactions = read.csv(args[2],sep=',',header=FALSE)
    } else {
        simulated_interactions = NULL
    }
} else {
    inFile <- 'data/Table-S1.csv'
    gsFile <- 'data/geneStats.csv'
    simulated_interactions = NULL
}
interactionFile <- 'data/interactions-present-in-data.txt'

inTab <- read.csv(inFile,stringsAsFactors=FALSE)
nr = nrow(inTab)

pThresh = 0.05

# Quantize the input data and write out the result.
# We will use Ckmeans to cluster each gene into nClusters
# ranges.
#nClusters <- 3
#nClusters <- 0

# Require at least 3 levels, so that the linear range of
# a gene's expression doesn't get squished into a single
# quantized level. If the expression values are 1.0, 100.0, 200.0,
# we want that to correspond to three levels, not one.
nClusters <- c(3,9)
geneStats = data.frame(fbid=character(),gene=character(),level=numeric(),min=numeric(),mean=numeric(),max=numeric(),sd=numeric())

nReplicates = floor((ncol(inTab)-3)/20)
outTab = emptyRow(nReplicates)

for (rr in 1:nr) {
    gvec <- as.numeric(inTab[rr,-1:-3])
    fbid = inTab[rr,1]
    fbnm = gname(fbid)
    ckms <- Ckmeans.1d.dp(gvec,nClusters)
    cFlags <- ckms['cluster'][[1]]

    # Compute the mean and sd within each level. We're going
    # to assume that within a level, expression values are
    # normally distributed.
    levels = unique(cFlags)
    if (length(levels) == 1) {
        cat("WTF? 1 level for ",fbid,"\n")
    }
    for(level in levels) {
        lVals = gvec[cFlags == level]
        sts = summary(lVals)
        sdev = sd(lVals)
        statRow = data.frame(fbid=fbid,gene=fbnm,level=level,min=sts['Min.'],mean=sts['Mean'],max=sts['Max.'],sd=sdev)
        geneStats = rbind(geneStats,statRow)
    }
    gName = inTab[rr,1]
    outTab[rr,1] = gName
    outTab[rr,-1:-3] = cFlags
}

#print("geneStats")
#print(geneStats)

colsPerCellType = nReplicates*5
colsPerRegion = nReplicates*4
colors.by.cellType = c(rep('blue',colsPerCellType),rep('green',colsPerCellType),rep('purple',colsPerCellType),rep('red',colsPerCellType))
colors.by.region = c(rep(c(rep('blue',nReplicates),rep('green',nReplicates),rep('black',nReplicates),rep('red',nReplicates),rep('orange',nReplicates)),4))
symbols.by.region = c(rep(c(rep('1',nReplicates),rep('2',nReplicates),rep('3',nReplicates),rep('4',nReplicates),rep('5',nReplicates)),4))

write.csv(outTab,paste(inFile,"-discrete-",nClusters[1],"-",nClusters[length(nClusters)],".csv",sep="",collapse=""))

row.names(geneStats) <- NULL
write.csv(geneStats,gsFile,row.names=FALSE)

interactions = read.table(interactionFile,header=FALSE,sep=",",stringsAsFactors=FALSE)

#i1 = interactions[1,]

test.interaction = function(i1,discreteDat,continDat) {
    #print("Testing interaction:")
    #print(i1)
    dRows = discreteDat[discreteDat[,1] == i1[1,1],]
    dRows = rbind(dRows,discreteDat[discreteDat[,1] == i1[1,2],])
    cRows = continDat[continDat[,1] == i1[1,1],-1:-3]
    cRows = rbind(cRows,continDat[continDat[,1] == i1[1,2],-1:-3])
    #print("cRows")
    #print(cRows)

    u1 = length(unique(as.numeric(dRows[1,4:length(dRows[1,])])))
    if (u1 == 1) {
        cat("SINGLE CLUSTER ",i1[1,1],": range",min(as.numeric(cRows[1,-1:-3])),"..",max(as.numeric(cRows[1,-1:-3])),"\n")
    }
    u2 = length(unique(as.numeric(dRows[2,4:length(dRows[2,])])))
    if (u2 == 1) {
        cat("SINGLE CLUSTER ",i1[1,2],": range",min(as.numeric(cRows[2,-1:-3])),"..",max(as.numeric(cRows[2,-1:-3])),"\n")
    }
    dRows2 = as.data.frame(t(dRows[,-1:-3]))
    conTab = table(dRows2)

    print("Contingency table:")
    print(conTab)
    cat("sum=",sum(conTab),"\n")

    if (all(dim(conTab) == 1)) {
        print("Aborting: 1x1 contingency table")
        return(1.0)
    }

    #X2 = chisq.test(conTab)
    #cat("ChiSq p-value",X2$p.value,"\n")

    FX2 = fun.chisq.test(conTab)
    #cat("FunChiSq p-value",FX2$p.value,"\n")

    g1 = i1[1,1]
    g2 = i1[1,2]
    #cat("p=",FX2$p.value,"\n")
    if (FX2$p.value <= pThresh) {
    #if (X2$p.value <= pThresh) {
        plot(as.numeric(cRows[1,-1:-3]),as.numeric(cRows[2,-1:-3]),col=colors.by.cellType,pch=symbols.by.region,xlab=gname(g1),ylab=gname(g2))
        tStr = paste(c(g1,gname(g1)," vs ",g2,gname(g2),": p=",signif(FX2$p.value,digits=5),"(T)"),sep=" ",collapse=" ")
        title(tStr)
        legend("topleft",legend=c("EC","EE","EB","ISC"),fill=c('blue','green','purple','red'))
        #plot(as.numeric(cRows[1,4:43]),as.numeric(cRows[2,4:43]),col=colors.by.region)
        #tStr = paste(c(g1," : ",g2,"p=",signif(FX2$p.value,digits=5),"(R)"),sep=" ",collapse=" ")
        #title(tStr)
        cat("Conserved interaction: ",tStr,"\n")
        print("Contingency table:")
        print(conTab)
        cat("\n")
    }

    return(FX2$p.value)
}

p.values = NULL
for (ii in 1:nrow(interactions)) {
#for (ii in 1:5) {
  interRow = interactions[ii,]
  p.values = c(p.values,test.interaction(interRow,outTab,inTab))
}
cat("p.values",p.values)

# Abandoned in favor of a manual BH correction.
#p.values = p.adjust(p.values,method="fdr")
#cat("adjusted p.values",p.values)

interactions[['p.value']] <- p.values
interactions[['BH.significant']] <- rep(FALSE,nrow(interactions))
interactions[['simulated']] <- rep(FALSE,nrow(interactions))

BH.alpha = 0.25

if (!is.null(simulated_interactions)) {
    for (jj in 1:nrow(interactions)) {
        for (kk in 1:nrow(simulated_interactions)) {
            if ((interactions[jj,1] == simulated_interactions[kk,1]) &&
                (interactions[jj,2] == simulated_interactions[kk,2])) {
                interactions[jj,'simulated'] <- TRUE
            }
        }
    }
}
sorted.interactions = interactions[order(interactions[,'p.value']),]

for (jj in 1:nrow(sorted.interactions)) {
    cat(jj,BH.alpha,nrow(sorted.interactions),sorted.interactions[jj,'p.value'],'<?',(jj*BH.alpha/nrow(sorted.interactions)),'\n')
    if (sorted.interactions[jj,'p.value'] < (jj*BH.alpha/nrow(sorted.interactions))) {
        sorted.interactions[jj,'BH.significant'] <- TRUE
    }
}

sorted.interactions[['n']] = 1:nrow(sorted.interactions)

print("Ajusted interaction significance:")
print(sorted.interactions)

