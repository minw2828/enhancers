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
infile1  <- args[index+1]  # 
infile2  <- args[index+2]  # 
outfile1 <- args[index+3]  # 
outfile2 <- args[index+4]  # 


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
# <- as.numeric(f22(outfile1, '[.]', 5)[1])
# <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
getSnp1Minus2 <- function(data1, data2) {
    tmp <- setdiff(data1$snpName, data2$snpName)
    return(res)
}

getResult <- function(data, tSep, outfile1, outfile2) {
    res1 <- getSnp1Minus2(data, histone1, histone2, tSep)
    write_file(res1, outfile1, FALSE)
    res2 <- getSnp1Minus2(data, histone2, histone1, tSep)
    write_file(res2, outfile2, FALSE)
}

run <- function() {
    data1 <- read_file(infile1)  
    data2 <- read_file(infile2)

    res1 <- setdiff(data1$snpName, data2$snpName)
    write_file(res1, outfile1, FALSE)

    res2 <- setdiff(data2$snpName, data1$snpName)
    write_file(res2, outfile2, FALSE)
}

run()



