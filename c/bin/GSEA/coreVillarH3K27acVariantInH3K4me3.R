# /usr/bin/rscript

## Description:
## This module checks if core H3K27ac variants are within H3K4me3 regions.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 19 December 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
#  <- as.numeric(f22(outfile1, '[.]', 5)[1])
# <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    dat <- fread(infile, header = FALSE)
    setnames(dat, cns)
    return(dat)
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

getBoolHistone <- function(data1, data2) {
    data1[, histone := 'NA']
    data1[snpName %in% data2$snpName, histone := unique(data2$histone)]
    return(data1[])
}

getRatioHistone <- function(data) {
    dt1 <- setDT(ddply(data, .(histone, pheno), summarize, countHistone = length(snpName)))
    dt2 <- setDT(ddply(data, .(pheno), summarize, countAll = length(snpName)))
    setkey(dt1, pheno)
    setkey(dt2, pheno)
    res <- merge(dt1, dt2, all.x = TRUE)
    res[, ratio := countHistone/countAll]
    res[, answer := paste(round(ratio * 100, 2), '%', sep = '')]
    return(res[])
}

write_file <- function(content, outfile) {
    Sep <- detechSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)

    data <- getBoolHistone(data1, data2)
    
    res <- getRatioHistone(data)

    write_file(res, outfile1)
}

run()



