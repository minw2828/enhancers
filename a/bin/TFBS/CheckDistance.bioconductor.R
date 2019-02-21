# /usr/bin/rscript

## Description:
## This script plots up a heatmap to show euclidean distance of dbSUPER 
## hits in the bovine genome. The mean position of a dbSUPER was used.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 29 March 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]
infile2  <- args[index+2]
outfile1 <- args[index+3]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables # 
round2 <- 2
#as.numeric(f22(outfile1, '[.]', 6)[1])
# <- f22(outfile1, '[.]', 4)[1]


# analysis #
rawProcess1 <- function(data1) {
    data1[, `:=` (newStart = ifelse(start <= end, start, end),
                  newEnd   = ifelse(start <= end, end, start))]
    data1[, `:=` (start = NULL, end = NULL)]
    setnames(data1, c('newStart', 'newEnd'), c('start', 'end'))
    return(data1[])
}

getNbMergePerElen <- function(data, elen, round2) {
    gr <- getGranges(data)
    grl <- split(gr, gr$db)
    rgrl <- reduce(grl, min.gapwidth = elen)
    nr <- elementLengths(grl)
    nrl <-  elementLengths(rgrl)
    dt1 <- data.table(db = names(nr), nrv = nr)
    dt2 <- data.table(db = names(nrl), nrlv = nrl)
    setkey(dt1, db)
    setkey(dt2, db)
    mdt <- merge(dt1, dt2)
    mdt[, `:=`(mergeDis = elen, 
               ratio2   = paste(round((as.numeric(nrv) - as.numeric(nrlv)) / as.numeric(nrv) * 100, round2), '%', sep = ''))]
    return(mdt[])
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)

    data <- rawProcess1(data1)

    res <- rbindlist(lapply(data2$elen, function(x) getNbMergePerElen(data, x, round2)))   
    write_file(res, outfile1, TRUE)
}


run()


