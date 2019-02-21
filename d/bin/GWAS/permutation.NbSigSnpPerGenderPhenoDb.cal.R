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
infile2  <- args[2]
outfile1 <- args[3]
# infile1  <- '/tmp//gender_pheno_snpName_effect_pvalue.10e-08.2016-01-19.a.tab'
# infile2  <- '/tmp//snpName_db.2016-01-19.a.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/permutation.NbSigSnpPerGenderPhenoDb/permutation.NbSigSnpPerGenderPhenoDb.10e-08.10000.6.2016-01-19.a.csv'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables # 
pSigLevel <- as.numeric(f22(outfile1, '[.]', 5)[1])
ntimes    <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

calNbSigSnp <- function(data, sigThreshold) {
   dat <- data[pvalue <= sigThreshold]
   res <- nrow(dat)
   return(res)
}

calNbSigSnpsDb <- function(data1, data2, sigThreshold) {
    setkey(data1, snpName)
    setkey(data2, snpName)
    mdt <- merge(data1, data2)
    res <- calNbSigSnp(mdt, sigThreshold)
    return(res)
}

## randomly draw $ndraw number of rows from data and calculate the number of significant SNPs within ##
drawSnps_calNbSigSnp <- function(data, ndraw, sigThreshold) {
    dat <- data[sample(.N, ndraw)]
    calNbSigSnp(dat, sigThreshold)
}

## calculate permutations ##
### n: The number of times of permutations ###
calculate_permutations <- function(data, ndraw, sigThreshold, n) {
    res <- replicate(n, drawSnps_calNbSigSnp(data, ndraw, sigThreshold), simplify = "vector")
    return(res)
}

getLoopOverList <- function(data1, data2) {
    genders <- unique(data1$gender)
    phenos  <- unique(data1$pheno)
    dbs     <- unique(data2$db)
    dt <- setDT(expand.grid(genders, phenos, dbs))
    setnames(dt, c('gender', 'pheno', 'db'))
    return(dt[])
}

get1result <- function(data1, data2, loopList, index, pSigLevel, ntimes) {
    oneRow <- loopList[index]
    dat1 <- data1[gender == oneRow$gender & pheno == oneRow$pheno]
    dat2 <- data2[db     == oneRow$db]
    ndraw <- nrow(dat2)
    NbSigSnpOriginal    <- calNbSigSnpsDb(dat1, dat2, pSigLevel)
    NbSigSnpsPermutaion <- calculate_permutations(dat1, ndraw, pSigLevel, ntimes)
    res <- data.table(NbSigSnp = c(NbSigSnpOriginal, NbSigSnpsPermutaion),
                      tpe      = c('ori', paste('permutation', seq(length(NbSigSnpsPermutaion)), sep = '')))
    res[, `:=` (gender = oneRow$gender,
                pheno  = oneRow$pheno,
                db     = oneRow$db)]
    return(res[])
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ",",
                row.names = FALSE, col.names = FALSE)
}

getResult <- function(data1, data2, loopList, pSigLevel, ntimes, outfile) {
    n <- nrow(loopList)
    res <- rbindlist(lapply(seq(n), function(x) get1result(data1, data2, loopList, x, pSigLevel, ntimes)))
    write_file(res, outfile)
}

run <- function() {
    data1 <- read_file(infile1) # all GWAS results 
    data2 <- read_file(infile2) # enhancer variants
    loopList <- getLoopOverList(data1, data2)
    getResult(data1, data2, loopList, pSigLevel, ntimes, outfile1)
}

run()



