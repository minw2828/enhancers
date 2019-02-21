# /usr/bin/rscript

## Description:
## This module plots up results from GSEA based on meta-analysis result.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 19 July 2016
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name> <infile1> <infile2> <infile3> <infile4> <outfile1> <outfile2> <outfile3> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1] # lesESnp results (H3K4me3 & H3K27ac)
outfile1 <- args[2] # csv


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pSigLvl <- as.numeric(f22(outfile1, '[.]', 4)[1])


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

getSigRatio <- function(data, pSigLvl) {
    dt1 <- setDT(ddply(data, .(pheno, db), summarize, nLesESnp = length(lesESnpName)))
    dt2 <- setDT(ddply(data[pmeta <= pSigLvl], .(pheno, db), summarize, nSigLesESnp = length(lesESnpName)))
    setkey(dt1, pheno, db)
    setkey(dt2, pheno, db)
    mdt <- merge(dt1, dt2, allow.cartesian = TRUE)
    mdt[, ratio := nSigLesESnp/nLesESnp]
    mdt[, answer := paste(round(ratio, 4) * 100, '%', sep='')]
    smdt <- mdt[order(db, pheno)]
    return(smdt)
}

write_file <- function(content, outfile) {
    Sep <- detechSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data <- read_file(infile1) 
    res <- getSigRatio(data, pSigLvl)
    write_file(res, outfile1)
}

run()



