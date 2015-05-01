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
    if (any(is.na(nm))) {
        nm = fbid
    }
    return(nm)
}


