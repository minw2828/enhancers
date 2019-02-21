# /usr/bin/rscript

## Description:
## This script liftOvers 1 source of mammalian enhancer regions to the bovine
## genome, and merges any output genomic intervals that were in <1bp to each 
## other, using Bioconductor.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 29 March 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <ofnpattern>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]  # 
outfile1 <- args[index+2]  # 


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
round2 <- 0


# analysis #
rawProcess1 <- function(data1) {
    data1[, `:=` (newStart = ifelse(start <= end, start, end),
                  newEnd   = ifelse(start <= end, end, start))]
    data1[, `:=` (start = NULL, end = NULL)]
    setnames(data1, c('newStart', 'newEnd'), c('start', 'end'))
    return(data1)
}

run <- function() {
    data1 <- read_file(infile1)

    data <- rawProcess1(data1)

    igr <- as(data, "GRanges")
    igrl <- split(igr, igr$db)
    rigrl <- reduce(igrl)
    
    dt <- cbind(data.table(db = names(elementLengths(rigrl)),
                           nb = elementLengths(rigrl),
                           s  = round(sd(width(rigrl)), round2)), 
                as.data.table(do.call(rbind, lapply(width(rigrl), function(x) do.call(rbind, list(round(summary(x), round2)))))))
    write_file(dt, outfile1, TRUE)
}

run()


