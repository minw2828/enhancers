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
outfile1 <- args[2]


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


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

filterByPSigLevel <- function(data, pSigLevel) {
    res <- data[pmeta <= pSigLevel]
    return(res)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = '\t',
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data <- read_file(infile1)
    res <- filterByPSigLevel(data, pSigLevel)
    write_file(res, outfile1)
}

run()


