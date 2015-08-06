# Compute input to CPX2 for detecting unique interactions.
#
# A single command-line argument is expected, which specifies
# a pattern (regular expression) to look for in column names.
# All columns that match the pattern will be pooled together;
# all the columns that do not match the pattern will pooled
# together in a second pool.
#
# Each two-gene interaction present in the Dutta dataset
# will be analyzed.
#
# For each interaction, the two pools described above are
# built. Trajectory files suitable for input to CPX2 are
# then written for each pool. The trajectory files are
# written to a subdirectory of ./work having the same
# name as the pattern given on the command line. The
# file names will be of the form
#
#  gene1-gene2-pattern.trj for the pattern pool;
#  gene1-gene2-__X__pattern.trj for the anti-pattern pool.
#
# A separate script runs CPX2 on the trajectory files
# and analyzes the output.

if (!exists("util.util",mode="function")) source("r/grn_util.r")

args <- commandArgs(trailingOnly=TRUE)
print(args)

interactionFile = "data/interactions-present-in-data.txt"

if (length(args) > 1) {
    inFile = args[2]
} else {
    inFile = 'data/Table-S1-discrete-3.csv'
}

dat = read.csv(inFile,stringsAsFactors=FALSE)
interactions = read.table(interactionFile,header=FALSE,sep=",",stringsAsFactors=FALSE)

# The pattern (regexp) that determines which columns we
# will pool against all other columns.
classPattern = args[1]

# Mitigate the maddening bizarreness of R's string
# handling facilities.
appendToAString = function(s,obj) {
    result = paste(s,as.character(obj),sep="",collapse="")
    return(result)
}

# Build a TRAJECTORY_VER2 string describing the network
# trajectory among the given data rows. Rows are gene
# IDs, columns are experimental conditions, timesteps, etc.
buildTrajectory = function(rows) {
    result = "TRAJECTORY_VER2\n"
    result = appendToAString(result,"1 2 0\n")
    result = appendToAString(result,"3\t3\n")
    result = appendToAString(result,rows[1,1])
    result = appendToAString(result,"\t")
    result = appendToAString(result,rows[2,1])
    result = appendToAString(result,"\n")
    result = appendToAString(result,"\n")
    result = appendToAString(result,ncol(rows)-1)
    result = appendToAString(result,"\n")
    for (col in 2:ncol(rows)) {
        result = appendToAString(result,rows[1,col]-1)
        result = appendToAString(result,"\t")
        result = appendToAString(result,rows[2,col]-1)
        result = appendToAString(result,"\n")
    }
    return(result)
}

# Select the columns that match or do not match a
# particular regexp.
selectClass = function(df,fbgns,colPattern,complement=FALSE) {
    result1 = df[df[,2] == fbgns[1],c(2,grep(colPattern,names(df),invert=complement))]
    result2 = df[df[,2] == fbgns[2],c(2,grep(colPattern,names(df),invert=complement))]
    result = rbind(result1,result2)
    result = result[,names(result) != "FBgn.1"]
    result = result[,names(result) != "X"]
    return(result)
}

# Compute the filename to which a trajectory will be written.
trajectoryFname = function(fbgns,colPattern,complement=FALSE) {
    geneNames = c(gname(fbgns[1,1]),gname(fbgns[1,2]))
    cat("Gene names for ",fbgns[1,1],fbgns[1,2]," = ",geneNames[1],geneNames[2],"\n")
    fName = paste(geneNames,collapse="-",sep="-")
    patPart = colPattern
    if (complement) {
        patPart = paste("__X__",patPart,collapse="",sep="")
    }
    fName = paste(fName,patPart,collapse="-",sep="-")
    fName = paste(fName,".trj",collapse="",sep="")
}

# Root of the project tree - BAD! Don't hard-code this,
# use an environment variable.
projRoot = "/home/jk/JKSync/jk/BINF/5353/GRNs"

# Write a trajectory to a file.
writeTrajectory = function(tStr,dirName,tFname) {
    mainDir = paste(projRoot,"/work",collapse="",sep="")
    subDir = dirName
    if (file.exists(file.path(mainDir,subDir))) {
        setwd(file.path(mainDir,subDir))
    } else {
        dir.create(file.path(mainDir,subDir))
        setwd(file.path(mainDir,subDir))
    }

    fconn = file(tFname)
    write(tStr,fconn)
    close(fconn)

    setwd(file.path(mainDir,".."))
}

# Process all interactions.
for (ii in 1:nrow(interactions)) {
    vinter = as.character(interactions[ii,])
    cat(paste(c("Processing interaction ",vinter,"\n"),sep=" ",collapse=" "))
    c1 = selectClass(dat,vinter,classPattern)
    c2 = selectClass(dat,vinter,classPattern,complement=TRUE)
    t1 = buildTrajectory(c1)
    t2 = buildTrajectory(c2)
    tFname1 = trajectoryFname(interactions[ii,],classPattern,complement=FALSE)
    tFname2 = trajectoryFname(interactions[ii,],classPattern,complement=TRUE)
    writeTrajectory(t1,classPattern,tFname1)
    writeTrajectory(t2,classPattern,tFname2)
}

cat("Done with the loop")

