# /usr/bin/rscript

## Description:
## This script summarizes fold change as per BH requires.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 15 December 2016
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
outfile1 <- args[2]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/permutation.NbSigSnpPerPhenoDb/NbSigSnp_tpe_pheno_db.10e-08.10000.2016-01-19.a.csv'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/permutation.NbSigSnpPerPhenoDb/pheno_db_rank.10e-08.10000.2016-01-19.a.csv'


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

getLoopOverList <- function(data) {
    dt <- unique(data[, .(pheno, db)])
    return(dt[])
}

getFoldChange <- function(data, loopList, index, ntimes) {
    oneRow <- loopList[index]
    dat <- data[pheno == oneRow$pheno & db == oneRow$db]
    sdt1 <- dat[tpe == 'ori']
    sdt2 <- dat[tpe != 'ori']
    fc <- sdt1$NbSigSnp / mean(sdt2$NbSigSnp)
    res <- data.table(pheno = oneRow$pheno, db = oneRow$db, foldChange = fc)
    return(res)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ',',
                row.names = FALSE, col.names = TRUE)
}

getResult <- function(data, loopList, ntimes, outfile) {
    n <- nrow(loopList)
    res <- rbindlist(lapply(seq(n), function(x) getFoldChange(data, loopList, x, ntimes)))
    write_file(res, outfile)
}

## run ##
run <- function() {
    data <- read_file(infile1) 
    loopList <- getLoopOverList(data)
    getResult(data, loopList, ntimes, outfile1)
}

run()



