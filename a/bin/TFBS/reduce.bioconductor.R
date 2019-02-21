# /usr/bin/rscript

## Description:
## This script returns non-overlapping genomic intervals and 
## 1000 Bull Genome Run 4 SNPs inside, using Bioconductor.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 07 April 2017
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf> <infile> <outfile>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]  # TFBS 
infile2  <- args[index+2]  # 1000 Bull Genome Run 4 SNPs 
outfile1 <- args[index+3]  # non-overlapping TFBS
outfile2 <- args[index+4]  # TFBS SNPs

# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
options(scipen = 999) # force not to use scienteific notation

# get global variables #
ORefseq  <- f22(outfile1, '[.]', 4)[1]


# analysis #
rawProcess1 <- function(data) {
    data[, `:=` (Start = ifelse(start <= end, start, end),
                 End   = ifelse(start <= end, end, start))]
    data[, `:=` (start = NULL, end = NULL)]
    setnames(data, c('Start', 'End'), c('start', 'end'))
    return(data)
}

getResult1 <- function(data) {
    rawProcess1(data)
    gr <- getGranges(data)
    rgr <- reduce(gr)
    tmp <- as.data.table(as.data.frame(rgr))
    res <- tmp[, .(seqnames, start, end)]
    return(res)
}

getResult2 <- function(data1, data2) {
    gr1 <- getGranges(data1)
    gr2 <- getGranges(data2)
    rgr1 <- reduce(gr1)
    tmp1 <- subsetByOverlaps(gr2, rgr1)
    tmp2 <- as.data.table(as.data.frame(tmp1))
    res <- tmp2[, .(snpName)]
    return(res)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)

    res1 <- getResult1(data1)
    res2 <- getResult2(data1, data2)

    write_file(res1, outfile1, FALSE)
    write_file(res2, outfile2, FALSE)
}


run()



