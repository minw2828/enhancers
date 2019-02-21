# /usr/bin/rscript

## Description:
## This module implements Gene Set Enrichment Analysis (GSEA) following step 1
## and step 2 in described reference below. Step 3 is unnecessary in our study.
##
## Reference:
## 2005, Subramanian A. et al., Gene set enrichment analysis: A knowledge-based
## approach for interpreting genome-wide expression profiles, PNAS
## 
## Interchange ideas between paper above and my study:
## Gene set: putatative promoter regions 
## ranking list: sorted minusLog10Pmeta
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
## Rscript <module_name> <infile1> <infile2> <outfile1> <outfile2> <outfile3> <outfile4>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1] # meta-analysis results 
infile2  <- args[2] # promoter snps 
outfile1 <- args[3] # GSEA results 
outfile2 <- args[4] # ES
outfile3 <- args[5] # lesESnp
outfile4 <- args[6] # esNull
outfile5 <- args[7] # EsRank


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
ntimes <- as.integer(f22(outfile2, '[.]', 4)[1])


# set seeds
set.seed(1234)


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getMinPmeta <- function(data) {
    data[, secondMinPmeta := tail(head(sort(unique(pmeta)), 2), 1)]
    data[, minPmeta := secondMinPmeta * 10e-08]
    data[pmeta == 0, pmeta := minPmeta]
    data[ ,`:=`(secondMinPmeta = NULL, minPmeta = NULL)]
    return(data)
}

getMinusLog10Pmeta <- function(data) {
    data[, minusLog10Pmeta := -log10(pmeta)]
    return(data)
}

mergeDatas <- function(data1, data2) {
    data1[snpName %in% data2$snpName, boolpromoter := TRUE]       # get promoter label
    data1[!(snpName %in% data2$snpName), boolpromoter := FALSE]   # get promoter label
    return(data1)
}

## calculate enrichment score ##
calculate_enrichment_score <- function(dat) {
    nFALSE <- nrow(dat[boolpromoter == FALSE])
    nR     <- sum(dat[boolpromoter == TRUE]$minusLog10Pmeta)
    dat[boolpromoter == TRUE, pHitMiss := minusLog10Pmeta/nR]
    dat[boolpromoter == FALSE, pHitMiss := -1/nFALSE]
    dat.sorted <- dat[order(minusLog10Pmeta, decreasing = TRUE)] # sort by -log10(pmeta), largest -> smallest
    dat.sorted[, accSum := cumsum(pHitMiss)]
    return(dat.sorted)
}

get_ES <- function(data) {
    res <- max(abs(data$accSum))
    if (res %in% data$accSum) {
        ES <- res
    } else {
        ES <- -res
    }
    return(ES)
}

## calculate ES(NULL) ##
calculate_esNull <- function(dat) {
    dat$boolpromoter <- sample(dat$boolpromoter)
    dat.sorted <- calculate_enrichment_score(dat)
    res <- get_ES(dat.sorted)
    return(res)
}

## calculate permutations ##
### n: The number of times of permutations ###
calculate_permutations <- function(data, ntimes) {
    res <- as.data.table(replicate(ntimes, calculate_esNull(data), simplify = "vector"))
    setnames(res, 'esNull')
    return(res)
}

## get leading edge subset ##
getLeadingEdgeSubset <- function(data) {
    ES <- get_ES(data)
    rowIndex <- data[,.I[data$accSum == ES]]
    res <- data[data[ , .I <= rowIndex] & boolpromoter == TRUE]
    return(res)
}

cal_ES_significance <- function(esNulls, ES) {
    esNulls[, id := paste('permutation', esNulls[,.I], sep = ":")]
    es <- data.table(esNull = ES, id = 'ES')
    mdt <- rbindlist(list(esNulls, es))
    mdt.sorted <- mdt[order(esNull)] # sort by esNull + ES
    nr.mdt.sorted <- nrow(mdt.sorted)
    ES.rank       <- mdt.sorted[,.I[mdt.sorted$id == 'ES']]
    if (ES.rank == nr.mdt.sorted) {
        res <- paste("<", round(1/nr.mdt.sorted, 6), sep = "")
    } else {
        res <- round((nr.mdt.sorted - ES.rank)/nr.mdt.sorted, 6)
    }
    return(res)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    getMinPmeta(data1)
    getMinusLog10Pmeta(data1)  
    data <- mergeDatas(data1, data2)
    data.es <- calculate_enrichment_score(data)
    write_file(data.es, outfile1)
    ES <- get_ES(data.es)
    write_file(ES, outfile2)
    lesESnp <- getLeadingEdgeSubset(data.es)
    write_file(lesESnp, outfile3)
    esNulls <- calculate_permutations(data.es, ntimes)
    write_file(esNulls, outfile4)  # must put here because later analysis will add columns
    ES.rank <- cal_ES_significance(esNulls, ES)
    write_file(ES.rank, outfile5)
}

run()



