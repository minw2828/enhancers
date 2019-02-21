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
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]
outfile3 <- args[5]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variable #
chrN   <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file1 <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

read_file2 <- function(infile) {
    tmp <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    cns <- c('phenoId', 'IntId', tmp[4:length(tmp)])
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

filterCrossAllCohorts <- function(data) {
    # filtered out is.na(effect), pvalue == 1 and pvalue == 0 
    tmp1 <- data[, if(all(!is.na(effect))) .SD, .(pheno, snpName)]
    tmp2 <- tmp1[, if(isTRUE(all(pvalue != 0))) .SD, .(pheno, snpName)]
    res  <- tmp2[, if(all(pvalue != 1)) .SD, .(pheno, snpName)]
    return(res)
}

get_z <- function(data) {
    data[, z := qnorm(pvalue)]
    return(data)
}

getStatistic_i <- function(data) {
    get_z(data)
    data[, se := abs(effect/z)]
    data[, weight := 1/se]
    data[, weighteffect := effect * weight]
    return(data)
}

assignGender <- function(data) {
    # The third ('sex') is 1 for cows and 2 for bulls.
    # The next three are for Fat, Milk and Prot Yield in that order.
    data[sex == 1, gender := 'cow']
    data[sex == 2, gender := 'bull']
    dat <- data[, .(gender, FY, MY, PY)]
    melted <- melt(dat, id.vars = "gender")
    setnames(melted, 'variable', 'pheno')
    return(melted)
}

getPhenoSummaries <- function(data) {
    data.summary <- setDT(ddply(data, .(gender, pheno), summarize, meanPheno = mean(value),
                                                                   sdPheno   = sd(value),
                                                                   minPheno  = min(value),
                                                                   maxPheno  = max(value),
                                                                   nAnimal   = length(value)))
    return(data.summary)
}

combineGenderPheno <- function(data1, data2) {
    # data1: meta-anlaysis results
    # data2: phenotype records
    setkey(data1, gender,pheno)
    setkey(data2, gender,pheno)
    mdt <- merge(data1, data2, all.x = TRUE)
    return(mdt)
}

addInfo2Data <- function(data1, data2) {
    data1.melted  <- assignGender(data1)
    data1.summary <- getPhenoSummaries(data1.melted)
    data <- combineGenderPheno(data2, data1.summary)
    return(data)
}

metaGWAS <- function(data) {
    dt <- setDT(aggregate(cbind(weighteffect, weight) ~ snpName + pheno, data, sum))
    dt[, snpeffectmeta := weighteffect / weight]
    dt[, npops := length(unique(data$gender))]
    dt[, varmeta := sqrt(npops/(weight * weight))]
    return(dt)
}

get_pvaluemeta <- function(data) {
    data[, zmeta := snpeffectmeta/varmeta]
    data[, pmeta := 2 * (1 - pnorm(abs(zmeta)))]
    return(data)
}

getMinPmeta <- function(data) {
    data[, secondMinPmeta := tail(head(sort(unique(pmeta)), 2), 1), by = .(pheno)]
    data[, minPmeta := secondMinPmeta * 10e-08]
    data[pmeta == 0, pmeta := minPmeta]
    data[ ,`:=`(secondMinPmeta = NULL, minPmeta = NULL)]
    return(data)
}

get_meta <- function(data) {
    res <- metaGWAS(data)
    get_pvaluemeta(res)
    getMinPmeta(res) # replace pmeta == 0 with second lowest pmeta/10e-08
    return(res[])
}

selectByFdr <- function(data, fdr) {
    dt <- setDT(ddply(data, .(pheno), summarize, fdrThreshold = compute.FDR(pmeta, fdr)))
    setkey(data, pheno)
    setkey(dt, pheno)
    mdt <- merge(data, dt, all.x = TRUE)
    res <- mdt[pmeta <= fdrThreshold]
    return(res)
}

get_nbSnps <- function(data, boolGender) {
    if (boolGender == TRUE) {
        dt <- setDT(aggregate(snpName ~ pheno + gender, data, length))
    } else {
        dt <- setDT(aggregate(snpName ~ pheno, data, length))
        dt[, gender := NA]
    }
    setnames(dt, 'snpName', 'nbSnps')
    res <- dt[, .(gender, pheno, nbSnps)]
    return(res[])
}

nbSnpLeft <- function(data1, data) {
    nbSnp_data1 <- get_nbSnps(data1, TRUE)
    nbSnp_data  <- get_nbSnps(data, TRUE)
    nbSnp_data1[, source := 'input']
    nbSnp_data[, source := 'getStatistics_i']
    res <- rbindlist(list(nbSnp_data1, nbSnp_data))
    return(res)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = '\t',
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data1 <- read_file1(infile1)
    data2 <- read_file2(infile2)
    data1.filtered <- filterCrossAllCohorts(data1)
    getStatistic_i(data1.filtered)
    data <- addInfo2Data(data2, data1.filtered)
    write_file(data, outfile1)
    meta <- get_meta(data)
    write_file(meta, outfile2)
    res <- nbSnpLeft(data1, data)
    write_file(res, outfile3)
}

run()


