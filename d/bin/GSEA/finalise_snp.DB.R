# /usr/bin/rscript

## Description:
## This script prepares bovine putative promoter regions to BH, request on 18 Jan 2016.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 11 March 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date>


# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
Db   <- f22(outfile1, '[.]', 6)[1]
Chew <- f22(outfile1, '[.]', 5)[1]


# analysis 
read_file1 <- function(infile, cols) {
    data <- fread(infile, header = FALSE, select = cols)
    cns  <- c("chr", "pos")
    setnames(data, cns)
    return(data)
}

read_file2 <- function(infile) {
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

## filter snps within putative promoters ##
select_by_chr <- function(data1, data2, chr_num) {
    dat1 <- data1[which(chr == chr_num),]
    dat2 <- data2[which(chr == chr_num),]
    promoter_positions <- unique(unlist(mapply(function(a, b) seq(a, b), dat2$start, dat2$end)))
    res <- dat1[dat1$pos %in% promoter_positions]
    return(res)
}

getResult <- function(data1, data2, Db, Chew) {
    result <- rbindlist(lapply(seq(30), function(x) select_by_chr(data1, data2, x)))
    result[, `:=` (db = Db, chew = Chew)]
    return(result[])
}

write_file <- function(content, outfile) {
    Sep <- detechSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data1 <- read_file1(infile1, c(3, 4))
    data2 <- read_file2(infile2)
    res <- getResult(data1, data2, Db, Chew)
    write_file(res, outfile1)
}

run()

