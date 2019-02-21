# /usr/bin/rscript

## Description:
## This module performs permutation test on meta-analysis results. 
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 14 July 2016
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]    # meta
infile2  <- args[index+2]    # dbSUPER merged SNPs
outfile1 <- args[index+3]
outfile2 <- args[index+4]


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
pheno     <- f22(outfile1, '[.]', 8)[1]
db        <- f22(outfile1, '[.]', 7)[1]
elen      <- as.numeric(f22(outfile1, '[.]', 6)[1])
pSigLevel <- as.numeric(f22(outfile1, '[.]', 5)[1])
ntimes    <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
calNbSigSnp <- function(data, sigThreshold) {
   dat <- data[pmeta <= sigThreshold]
   res <- nrow(dat)
   return(res)
}

calNbSigSnpsDb <- function(data1, data2, sigThreshold) {
    setkey(data1, snpName)
    setkey(data2, snpName)
    mdt <- merge(data1, data2)
    res <- calNbSigSnp(mdt, sigThreshold)
    return(res)
}

## randomly draw $ndraw number of rows from data and calculate the number of significant SNPs within ##
drawSnps_calNbSigSnp <- function(data, ndraw, sigThreshold) {
    dat <- data[sample(.N, ndraw, replace = TRUE)]
    calNbSigSnp(dat, sigThreshold)
}

## calculate permutations ##
### n: The number of times of permutations ###
calculate_permutations <- function(data, ndraw, sigThreshold, n) {
    res <- replicate(n, drawSnps_calNbSigSnp(data, ndraw, sigThreshold), simplify = "vector")
    return(res)
}

getFoldChange <- function(data, ntimes) {
    sdt1 <- data[tpe == 'ori']
    sdt2 <- data[tpe != 'ori']
    fc <- sdt1$NbSigSnp / mean(sdt2$NbSigSnp)
    return(fc)
}

getRanking <- function(data, ntimes) {
    data[, Rank := .I]
    data[, rank := ifelse(tpe == 'ori' & Rank == ntimes + 1, paste('<', 1/ntimes, sep = ''), as.character(Rank/ntimes))]
    res <- unique(data[tpe == 'ori', .(rank)])
    return(res)
}

getResult <- function(data1, data2, pSigLevel, ntimes, outfile1, outfile2) {
    ndraw <- nrow(data2)
    NbSigSnpOriginal    <- calNbSigSnpsDb(data1, data2, pSigLevel)
    NbSigSnpsPermutaion <- calculate_permutations(data1, ndraw, pSigLevel, ntimes)
    res1 <- data.table(NbSigSnp = c(NbSigSnpOriginal, NbSigSnpsPermutaion),
                       tpe      = c('ori', paste('permutation', seq(length(NbSigSnpsPermutaion)), sep = '')))
    write_file(res1, outfile1, FALSE)
    res2 <- data.table(rank       = getRanking(res1, ntimes),
                       foldChange = getFoldChange(res1, ntimes))
    write_file(res2, outfile2, TRUE)
}

run <- function() {
    data1 <- read_file(infile1) 
    data2 <- read_file(infile2) 

    getResult(data1, data2, pSigLevel, ntimes, outfile1, outfile2)
}

run()



