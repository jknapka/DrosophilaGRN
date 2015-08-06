args <- commandArgs(trailingOnly=TRUE)

print(args)
print(args[1])
fname <- args[1]
print(fname)

inDat <- read.csv(fname,stringsAsFactors=FALSE,header=FALSE)

gnameFile <- 'data/Dutta_FlyBase_IDs.txt'
gnameTab <- read.csv(gnameFile,sep='\t',stringsAsFactors=FALSE)

gname <- function(fbid) {
    nm = gnameTab[gnameTab[,2] == fbid,4]
    if (any(is.na(nm))) {
        nm = fbid
    }
    return(nm)
}

for (rr in 1:nrow(inDat)) {
  for (cc in 1:ncol(inDat)) {
      val = inDat[rr,cc]
      wtf = grepl(pattern='^FB',x=val)
      # cat(c("Looking at ",rr,cc,val,wtf,'\n'))
      if (wtf) {
          inDat[rr,cc] = gname(val)
      }
  }
}

write.table(inDat,file="",sep=",")

