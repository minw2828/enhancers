# /usr/bin/rscript

## Description:
## This script prepares bovine putative promoter regions to BH, request on 18 Jan 2016.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 04 March 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile> <outfile>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
outfile1 <- args[2]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
library(GenomicRanges); library(Biostrings); library(GenomicAlignments); library(GenomicFeatures)
library(data.table); library(plyr); library(dplyr)
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
Db   <- f22(outfile1, '[.]', 6)[1]
Chew <- f22(outfile1, '[.]', 5)[1]


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

detechSep <- function(file) {
    fmt <- f22(file, '[.]', 1)
    if (fmt == 'csv') {
        return(',')
    } else if (fmt %in% c('tab', 'bed', 'sam')) {
        return('\t')
    } else {
        return (' ')
    }
}

getNonduplicate <- function(data, Db, Chew) {
    gr  <- as(data, 'GRanges')
    rgr <- reduce(gr)
    dt  <- as.data.table(as(rgr, 'data.frame'))
    sdt <- dt[order(seqnames, start, end)]
    sdt[, `:=` (db = Db, chew = Chew)]
    return(sdt[, .(seqnames, start, end, db, chew)])
}

write_file <- function(content, outfile) {
    Sep <- detechSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data <- read_file(infile1)
    res <- getNonduplicate(data, Db, Chew)
    write_file(res, outfile1)
}

run()



