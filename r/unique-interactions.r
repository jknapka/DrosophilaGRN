args <- commandArgs(trailingOnly=TRUE)
print(args)

dat = read.csv('data/Table-S1-discrete-3.csv',stringsAsFactors=FALSE)

classPattern = args[1]

buildTrajectory = function(rows) {
    cat(paste(colnames(rows),collapse=",",sep=","))
    for (rr in 1:nrow(rows)) {
    }
}

selectClass = function(df,fbgns,colPattern,complement=FALSE) {
    result = df[df[,2]==fbgns,c(2,grep(colPattern,names(df),invert=complement))]
    return(result)
}


