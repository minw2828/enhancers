# /usr/bin/rscript

## Description:
## This script performs meta-analysis by re-engineering BH's zScore.f90 script.
##   
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 06 July 2016 
## 
## Date modified and reason: 
## 13 July 2016    Use qnorm(value) to calcualte x value.
##
## Execution: 
## Rscript <module_name> 


# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variable #
pSigLevel <- as.numeric(f22(outfile1, '[.]', 5)[1])
chrN      <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getData <- function(data1, data2) {
    setkey(data1, snpName)   
    setkey(data2, snpName)
    mdt <- merge(data1, data2, allow.cartesian = TRUE)
    return(mdt)
}

countSigSnp <- function(data, pSigLevel) {
    dat <- data[pmeta <= pSigLevel]
    dt <- setDT(ddply(dat, .(pheno, correction), summarize, nbSig = length(snpName)))
    res <- spread(dt, correction, nbSig)   
    return(res[, .(pheno, before, after)])
}

countSigESnp <- function(data, pSigLevel) {
    dat <- data[pmeta <= pSigLevel]
    dt <- setDT(ddply(dat, .(pheno, db, correction), summarize, nbSig = length(snpName)))
    res <- spread(dt, correction, nbSig)
    return(res[, .(pheno, db, before, after)])
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ',',
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)

    data <- getData(data1, data2)
   
    res1 <- countSigSnp(data1, pSigLevel)
    write_file(res1, outfile1)

    res2 <- countSigESnp(data, pSigLevel)
    write_file(res2, outfile2)
}

run()


