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
infile1  <- args[1]
outfile1 <- args[2]


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pSig <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    dat <- fread(infile, header = FALSE)
    setnames(dat, cns)
    return(dat)
}

countSigPerChrPheno <- function(data, pSig) {
    dat <- data[pmeta <= pSig]
    dt <- setDT(ddply(dat, .(chr, pheno), nrow))
    setnames(dt, 'V1', 'count')
    return(dt)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ",",
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data <- read_file(infile1)
    res <- countSigPerChrPheno(data, pSig)
    write_file(res, outfile1)
}

run()



