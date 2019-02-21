# /usr/bin/rscript

## Description:
## This script finds overlap between SNPs
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 25 November 2016
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name>  <infile1> <outfile1>


# get input arguments #
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
pmetaSig <- f22(outfile1, '[.]', 5)[1]
round2   <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getBoolenhancerPerDb <- function(data1, data2) {
    data1[snpName %in% data2$snpName, boolEnhancer := TRUE]
    data1[!(snpName %in% data2$snpName), boolEnhancer := FALSE]
    return(data1)
}

getSigSnpPerDb <- function(data1, data2, pmetaSig, Db) {
    pmetaSig <- as.numeric(pmetaSig)
    dat1 <- data1[pmeta <= pmetaSig]
    dat2 <- data2[db == Db]
    getBoolenhancerPerDb(dat1, dat2)
    res <- dat1[boolEnhancer == TRUE]
    res[, db := Db]
    return(res[])
}

transformData <- function(data1, data2, pmetaSig, dbs) {
    X <- rbindlist(lapply(dbs, function(x) getSigSnpPerDb(data1, data2, pmetaSig, x)))
    return(X)
}

getRatio <- function(data, DB1, DB2) {
    # overlaps(DB1, DB2)/DB1
    dat1 <- data[db == DB1]
    dat2 <- data[db == DB2]
    setkey(dat1, snpName,pheno)
    setkey(dat2, snpName,pheno)
    mdt <- merge(dat1, dat2, allow.cartesian = TRUE)
    res <- nrow(mdt)/nrow(dat1)
    dt <- data.table(db1 = DB1, db2 = DB2, ratio = res)
    return(dt)
}

getResult <- function(data, dbs, round2) {
    dt <- setDT(expand.grid(dbs, dbs))
    ns <- seq(nrow(dt))
    tmp <- rbindlist(lapply(ns, function(x) getRatio(data, dt[x]$Var1, dt[x]$Var2)))
    tmp[, ratio := paste(round(ratio, round2)*10**2, '%', sep = '')]
    res <- reshape(tmp, timevar = 'db2', idvar = 'db1', direction = 'wide')
    setnames(res, c('db', as.character(res$db)))
    setnames(res, 'db', 'NA')
    return(res)
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

write_file <- function(content, outfile) {
    Sep <- detechSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)

    dbs <- unique(data2$db)

    Data <- transformData(data1, data2, pmetaSig, dbs)

    res <- getResult(data, dbs, round2)
    write_file(res, outfile1)
}

run ()

