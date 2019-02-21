# /usr/bin/rscript

## Description:
## This module performs permutation test on meta-analysis results. 
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 30 March 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]    # dbSUPER region
infile2  <- args[index+2]    # 1000 Bull Genome SNPs
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
elen <- as.numeric(f22(outfile1, '[.]', 6)[1])
db   <- f22(outfile1, '[.]', 5)[1]


# analysis #
getSnpsDb <- function(gr1, gr2) {
    tmp1 <- subsetByOverlaps(gr2, gr1)
    tmp2 <- as.data.table(as.data.frame(tmp1))
    res <- tmp2[, .(snpName)]
    return(res)
}

run <- function() {
    data1 <- read_file(infile1) 
    data2 <- read_file(infile2) 

    gr1 <- reduce(getGranges(data1), min.gapwidth = elen)
    gr2 <- getGranges(data2)
    
    res <- getSnpsDb(gr1, gr2)
    write_file(res, outfile1, FALSE)
}

run()



