# /usr/bin/rscript

## Description:
## This module computes false discovery rate (FDR) for input datatable with target vector.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 08 September 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1] # lesESnp
infile2  <- args[2] # meta-analysis 
infile3  <- args[3] # enhancer SNPs
outfile1 <- args[4]


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    dat <- fread(infile, header = FALSE)
    setnames(dat, cns)
    return(dat)
}

getNb <- function(data1, data2, data3) {
    dt1 <- setDT(ddply(data1, .(pheno, db), summarize, nbLesESnp = length(lesESnpName)))
    dt2 <- setDT(ddply(data2, .(pheno), summarize, nbSnp = length(snpName)))
    dt3 <- setDT(ddply(data3, .(db), summarize, nbESnp = length(snpName)))
    dt2[, nbESnp := dt3$nbESnp]
    setkey(dt1, pheno)
    setkey(dt2, pheno)
    mdt <- merge(dt1, dt2)
    return(mdt)
} 

getRatio <- function(data) {
    data[, ratioLesESnpPerESnp := nbLesESnp/nbESnp]
    data[, ratioLesESnpPerSnp := nbLesESnp/nbSnp]
    data[, ratioESnpPerSnp := nbESnp/nbSnp]
    data[, answerLesESnpPerESnp := paste(round(ratioLesESnpPerESnp, 4)*100, '%', sep = '')]
    data[, answerLesESnpPerSnp := paste(round(ratioLesESnpPerSnp, 4)*100, '%', sep = '')]
    data[, answerESnpPerSnp := paste(round(ratioESnpPerSnp, 4)*100, '%', sep = '')]
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ",",
                row.names = FALSE, col.names = FALSE)
}


run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)    
    data3 <- read_file(infile3)

    data <- getNb(data1, data2, data3)

    getRatio(data)
    write_file(data, outfile1)
}

run()



