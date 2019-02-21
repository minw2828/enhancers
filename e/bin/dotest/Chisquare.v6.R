# /usr/bin/rscript

## Description:
## On 14 September 2017, MG, AC, TPH and MW had a skype call to BH. MG raises
## a concern on the permutation test that MW performed on her first paper.
## MG's concern was that whether the enrichment signal from H3K4me3 regions
## was due to LD to mutation in gene, or true enrichment within H3K4me3 regions.
## To address the concern, MG proposes the following analyses that will be
## performed in this folder.
##
## This script performs a Chi-Square test as MG suggests on 18 September.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 18 September 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <ofnpattern>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1   <- args[index+1]    # result of makeTable
outfile1  <- args[index+2]    # 
outfile2  <- args[index+3]    #


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
outpath  <- file.path(path_pre, 'analyses', date, psf, 'out')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
source(file.path(binpath, 'getData1234.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables #
cov1     <- as.numeric(f22(outfile1, '[.]', 6)[1])    # exclusive
cov2     <- as.numeric(f22(outfile1, '[.]', 5)[1])    # exclusive
ORefseq  <- f22(outfile1, '[.]', 4)[1]


# get supporting functions #
#Task <- f21(f22(program, '/', 1), '[.]', 1)


# set seed #
set.seed(1234)    # Need to set seed when using caret package


# analysis #
rawProcess <- function(data, cov1, cov2) {
    dat <- data[ratioOfProInDb > cov1 & ratioOfProInDb < cov2]
    dt <- dat[!(nbQltInDb == 0 & nbQltOutDb == 0)]
    setnames(dt, c('nbQltInDb', 'nbQltOutDb'), c('obsNbQtlInDb', 'obsNbQtlOutDb'))
    return(dt[])
} 

getCols <- function(data) {
    data[, totalNbSipSnp := obsNbQtlInDb + obsNbQtlOutDb]
    tmp <- colSums(data[, .(obsNbQtlInDb, obsNbQtlOutDb)])
    total <- sum(data$totalNbSipSnp)
    data[, `:=` (expNbQtlInDb  = tmp['obsNbQtlInDb']  / total * totalNbSipSnp,
                 expNbQtlOutDb = tmp['obsNbQtlOutDb'] / total * totalNbSipSnp)]
    return(data)
}

getResult <- function(data) {
    if (nrow(data) != 0) {
        data[, chi := (obsNbQtlInDb - expNbQtlInDb) ** 2 / expNbQtlInDb]
        chisqv <- sum(data$chi)
        df <- nrow(data) - 1
        chisqp <- pchisq(chisqv, df, lower.tail = FALSE)
        res <- data.table(chisqValue      = chisqv, 
                          chisqPvalue     = chisqp, 
                          degreeOfFreedom = df, 
                          tpe             = 'BrianCaffo.chisqtest')
    } else {
        res <- data.table(chisqValue      = NA, 
                          chisqPvalue     = NA,
                          degreeOfFreedom = NA,
                          tpe             = 'BrianCaffo.chisqtest')
    }
    return(res)
}

run <- function() {
    data1 <- read_file(infile1)

    data <- rawProcess(data1, cov1, cov2)
    getCols(data)
    write_file(data, outfile1, FALSE)

    res <- getResult(data)
    write_file(res, outfile2, FALSE)
}

run()





