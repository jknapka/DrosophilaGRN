# Utility functions.

util.util <- function() {
  # I exist only as a flag.
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

