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


# get global variables # 
# <- as.numeric(f22(outfile1, '[.]', 5)[1])
# <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

rawProcess <- function(data, tSep) {
    data[, snpName := paste(chr, pos, sep = tSep)]
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

getOverlapped <- function(data) {
    tmp <- Reduce(intersect, data[, .(list(unique(snpName))), db]$V1)
    dbs <- unique(data$db)
    res <- data.table(snpName = tmp, db = paste(unique(unlist(strsplit(dbs, ':'))), collapse = "_"))
    return(res)
}

write_file <- function(content, outfile) {
    Sep <- detechSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data <- read_file(infile1)  
    res <- getOverlapped(data)
    write_file(res, outfile1)
}

run()



