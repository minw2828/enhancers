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
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]  # meta
infile2  <- args[index+2]  # TFBS
outfile1 <- args[index+3]  
outfile2 <- args[index+4]
outfile3 <- args[index+5]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
options(scipen = 999) # force not to use scienteific notation

# get global variables #
round2    <- 2
pSigLevel <- as.numeric(f22(outfile1, '[.]', 4)[1])


# read file #
getData <- function(data1, data2) {
    setkey(data1, chr,pos,snpName)
    setkey(data2, chr,pos,snpName)
    mdt <- merge(data1, data2, all.y = TRUE)
    return(mdt)
}

countSigSnpPerPheno <- function(data, pSigLevel, round2) {
    dt1 <- setDT(ddply(data[pmeta <= pSigLevel], .(pheno), summarize, nbSigSnp = length(snpName)))
    dt2 <- setDT(ddply(data, .(pheno), summarize, nbSnp = length(snpName)))
    setkey(dt1, pheno)
    setkey(dt2, pheno)
    mdt <- merge(dt1, dt2)
    mdt[, ratio := nbSigSnp/nbSnp]
    mdt[, answer := paste(round(ratio * 100, round2), '%', sep = '')]
    return(mdt[])
}

countSigSnpPerPhenoChr <- function(data, pSigLevel) {
    dat <- data[pmeta <= pSigLevel]
    dt <- setDT(ddply(dat, .(pheno, chr), summarize, nbSnp = length(snpName)))
    dt.ordered <- dt[order(pheno, nbSnp, chr, decreasing = TRUE)]
    return(dt.ordered)
}

countSigSnpPerPhenoChrDb <- function(data1, data2, pSigLevel) {
    data <- getData(data1, data2)
    dat <- data[pmeta <= pSigLevel]
    dt <- setDT(ddply(dat, .(pheno, db, chr), summarize, nbSnp = length(snpName)))
    dt.ordered <- dt[order(pheno, db, nbSnp, chr)]
    return(dt.ordered)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)

    res1 <- countSigSnpPerPheno(data1, pSigLevel, round2)
    write_file(res1, outfile1, TRUE)

    res2 <- countSigSnpPerPhenoChr(data1, pSigLevel)
    write_file(res2, outfile2, TRUE)

    res3 <- countSigSnpPerPhenoChrDb(data1, data2, pSigLevel)
    write_file(res3, outfile3, TRUE)
}

run()


