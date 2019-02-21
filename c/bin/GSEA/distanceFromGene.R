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
infile1  <- args[1] # p-valule 
infile2  <- args[2] # gene region 
infile3  <- args[3] # enhancers 
outfile1 <- args[4]


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
chrN      <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getStats <- function(data1, data2, data3, pSigLevel, DB) {
    dat3 <- data3[db == DB]
    dat1 <- data1[pmeta <= pSigLevel & snpName %in% dat3$snpName]
    dat1[, `:=` (front = data2$start - pos,
                 back  = pos - data2$end)]
    res1 <- dat1[pos >= data2$start & pos <= data2$end]
    res2 <- dat1[front > 0][which.min(front)]
    res3 <- dat1[back > 0][which.min(back)]
    res4 <- dat1[front > 0][which.max(front)]
    res5 <- dat1[back > 0][which.max(back)]
    dt <- data.table(db             = rep(DB, 5),
                     snpName     = c(ifelse(nrow(res1) == 0, NA, res1$snpName), res2$snpName, res3$snpName, res4$snpName, res5$snpName),  
                     pheno          = c(ifelse(nrow(res1) == 0, NA, res1$pheno), res2$pheno, res3$pheno, res4$pheno, res5$pheno), 
                     snpeffectmeta  = c(ifelse(nrow(res1) == 0, NA, res1$snpeffectmeta), res2$snpeffectmeta, res3$snpeffectmeta, res4$snpeffectmeta, res5$snpeffectmeta), 
                     pmeta          = c(ifelse(nrow(res1) == 0, NA, res1$pmeta), res2$pmeta, res3$pmeta, res4$pmeta, res5$pmeta), 
                     minFront2Gene = c(ifelse(nrow(res1) == 0, NA, 0), res2$front, NA, NA, NA),
                     maxFront2Gene = c(ifelse(nrow(res1) == 0, NA, 0), NA, res4$front, NA, NA),
                     minBack2Gene  = c(ifelse(nrow(res1) == 0, NA, 0), NA, NA, res3$back, NA),
                     maxBack2Gene  = c(ifelse(nrow(res1) == 0, NA, 0), NA, NA, NA, res5$back))
    return(dt)
}

getResult <- function(data1, data2, data3, pSigLevel) {
    dbs <- unique(data3$db)
    res <- rbindlist(lapply(dbs, function(x) getStats(data1, data2, data3, pSigLevel, x)))
    return(res)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ',',
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    data3 <- read_file(infile3)

    res <- getResult(data1, data2, data3, pSigLevel)
    write_file(res, outfile1)
}

run()


