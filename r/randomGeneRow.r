if (!exists("util.util",mode="function")) {
    # print("Reading grn_util.r")
    source("r/grn_util.r")
}

geneStats = read.csv('data/geneStats.csv',sep=",",header=TRUE)

randomGeneRow = function(n) {
    result = emptyRow()

    # First, pick a gene.
    genes = levels(geneStats[,'gene'])
    gene = sample(genes,1)
    grows = geneStats[geneStats[,'gene']==gene,]

    geneNm = paste("FBid100",n,sep="",collapse="")
    cgNum = paste("CG10",n,sep="",collapse="")
    result[1,1] = geneNm
    result[1,2] = cgNum
    result[1,3] = geneNm
    
    glevs = grows[,'level']
    # We need 20 pairs of random samples.
    for (ss in 1:20) {
        # Pick a level.
        glev = sample(glevs,1)
        grow = grows[grows[,'level']==glev,]
        for (gss in 1:2) {
            mean = grow[,'mean']
            sdev = grow[,'sd']
            if (is.na(sdev) || is.nan(sdev)) {
                #print("CHANGING NaN sdev to 0.001")
                sdev = 0.001
            }
            expVal = rnorm(1,mean,sdev)
            if (is.nan(expVal)) {
                cat("NaN for ",as.character(grow)," mean=",mean," sd=",sdev,"\n")
            }
            result[1,2*(ss-1)+gss+3] = expVal
        }
    }

    return(result)
}

emptyRow = function() {
    nr = 0
    outTab <- data.frame(FBgn = rep("",nr),
                         GNum = rep("",nr),
                         GName = rep("",nr),
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
    colnames(outTab)[2] <- "CG number"
    colnames(outTab)[3] <- "Gene ID" 
    return(outTab)
}

generateRandomRows <- function(nGenes) {
    geneStats = read.csv('data/geneStats.csv',header=TRUE)

    nGenes = as.numeric(args[1])
    result = emptyRow()
    for (gn in 1:nGenes) {
        rg = randomGeneRow(gn)
        result = rbind(result,rg)
    }

    return(result)
}

fillDataTable <- function(outFile) {
    tab = read.csv('data/Table-S1.csv')
    for (ss in 1:nrow(tab)) {
        rr = randomGeneRow(0)
        tab[ss,-1:-3] = rr[1,-1:-3]
    }
    write.csv(tab,outFile,row.names=FALSE)
}


args <- commandArgs(trailingOnly=TRUE)
fillDataTable(args[1])

if (FALSE) {
    print(generateRandomRows(as.numeric(args[1])))
}

