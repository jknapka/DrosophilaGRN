require('Ckmeans.1d.dp')
require('FunChisq')

if (!exists("util.util",mode="function")) source("r/grn_util.r")

inFile <- 'data/Table-S1.csv'
interactionFile <- 'data/interactions-present-in-data.txt'

inTab <- read.csv(inFile,stringsAsFactors=FALSE)
nr = nrow(inTab)

pThresh = 0.05

# Quantize the input data and write out the result.
# We will use Ckmeans to cluster each gene into nClusters
# ranges.
nClusters <- 3
outTab <- data.frame(FBgn = rep("",nr),
                     EC_R1_1 = rep(NA,nr),
                     EC_R1_2 = rep(NA,nr),
                     EC_R2_1 = rep(NA,nr),
                     EC_R2_2 = rep(NA,nr),
                     EC_R3_1 = rep(NA,nr),
                     EC_R3_2 = rep(NA,nr),
                     EC_R4_1 = rep(NA,nr),
                     EC_R4_2 = rep(NA,nr),
                     EC_R5_1 = rep(NA,nr),
                     EC_R5_2 = rep(NA,nr),
                     EB_R1_1 = rep(NA,nr),
                     EB_R1_2 = rep(NA,nr),
                     EB_R2_1 = rep(NA,nr),
                     EB_R2_2 = rep(NA,nr),
                     EB_R3_1 = rep(NA,nr),
                     EB_R3_2 = rep(NA,nr),
                     EB_R4_1 = rep(NA,nr),
                     EB_R4_2 = rep(NA,nr),
                     EB_R5_1 = rep(NA,nr),
                     EB_R5_2 = rep(NA,nr),
                     EE_R1_1 = rep(NA,nr),
                     EE_R1_2 = rep(NA,nr),
                     EE_R2_1 = rep(NA,nr),
                     EE_R2_2 = rep(NA,nr),
                     EE_R3_1 = rep(NA,nr),
                     EE_R3_2 = rep(NA,nr),
                     EE_R4_1 = rep(NA,nr),
                     EE_R4_2 = rep(NA,nr),
                     EE_R5_1 = rep(NA,nr),
                     EE_R5_2 = rep(NA,nr),
                     ISC_R1_1 = rep(NA,nr),
                     ISC_R1_2 = rep(NA,nr),
                     ISC_R2_1 = rep(NA,nr),
                     ISC_R2_2 = rep(NA,nr),
                     ISC_R3_1 = rep(NA,nr),
                     ISC_R3_2 = rep(NA,nr),
                     ISC_R4_1 = rep(NA,nr),
                     ISC_R4_2 = rep(NA,nr),
                     ISC_R5_1 = rep(NA,nr),
                     ISC_R5_2 = rep(NA,nr),
                     stringsAsFactors=FALSE)
for (rr in 1:nr) {
    gvec <- as.numeric(inTab[rr,4:43])
    ckms <- Ckmeans.1d.dp(gvec,nClusters)
    cFlags <- ckms['cluster'][[1]]
    gName = inTab[rr,1]
    outTab[rr,1] = gName
    outTab[rr,2:41] = cFlags[1:40]
    #for (ss in 2:41) {
    #    cf = cFlags[ss-1]
    #    outTab[rr,ss] <- cf
    #}
}

colors.by.cellType = c(rep('blue',10),rep('green',10),rep('purple',10),rep('red',10))
colors.by.region = c(rep(c('blue','blue','green','green','black','black','red','red','orange','orange'),4))
symbols.by.region = c(rep(c('1','1','2','2','3','3','4','4','5','5'),4))

write.csv(outTab,paste("data/Table-S1-discrete-",nClusters,".csv",sep="",collapse=""))

interactions = read.table(interactionFile,header=FALSE,sep=",",stringsAsFactors=FALSE)

#i1 = interactions[1,]

test.interaction = function(i1,discreteDat,continDat) {
    #print("Testing interaction:")
    #print(i1)
    dRows = discreteDat[discreteDat[,1] == i1[1,1],]
    dRows = rbind(dRows,discreteDat[discreteDat[,1] == i1[1,2],])
    cRows = continDat[continDat[,1] == i1[1,1],]
    cRows = rbind(cRows,continDat[continDat[,1] == i1[1,2],])
    #print(dRows)

    conTab = matrix(nrow=nClusters,ncol=nClusters,0)
    for (cc in 2:ncol(dRows)) {
        count = conTab[dRows[1,cc],dRows[2,cc]]
        count = count+1
        conTab[dRows[1,cc],dRows[2,cc]]=count
    }

    #print("Contigency table:")
    #print(conTab)

    X2 = chisq.test(conTab)
    #cat("ChiSq p-value",X2$p.value,"\n")

    FX2 = fun.chisq.test(conTab)
    #cat("FunChiSq p-value",FX2$p.value,"\n")

    g1 = i1[1,1]
    g2 = i1[1,2]
    if (FX2$p.value <= pThresh) {
        plot(as.numeric(cRows[1,4:43]),as.numeric(cRows[2,4:43]),col=colors.by.cellType,pch=symbols.by.region,xlab=gname(g1),ylab=gname(g2))
        tStr = paste(c(gname(g1)," vs ",gname(g2),": p=",signif(FX2$p.value,digits=5),"(T)"),sep=" ",collapse=" ")
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

for (ii in 1:nrow(interactions)) {
#for (ii in 1:5) {
  interRow = interactions[ii,]
  test.interaction(interRow,outTab,inTab)
}

