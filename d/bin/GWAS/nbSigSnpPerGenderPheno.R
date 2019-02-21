# /usr/bin/rscript

## Description:
## This script calculates the number of SNP per gender, phenotype cohort.
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
## Rscript <module_name> <infile1> <outfile1> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]  # GWAS result 
outfile1 <- args[2]
# infile1  <- '/tmp//gender_pheno_chr_pos_effect_pvalue.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/Manhattan.annotate_snp_gene_enhancer.NACohorts/output.10.2.a.2016-01-19.a.png'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pSigLevel <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

countSigSnpPerGenderPheno <- function(data, pSigLevel, round2) {
    dt1 <- setDT(ddply(data[pvalue <= pSigLevel], .(gender, pheno), summarize, nbSigSnp = length(snpName)))
    dt2 <- setDT(ddply(data, .(gender, pheno), summarize, nbSnp = length(snpName)))
    setkey(dt1, gender,pheno)
    setkey(dt2, gender,pheno)
    mdt <- merge(dt1, dt2)
    mdt[, ratio := nbSigSnp/nbSnp]
    mdt[, answer := paste(round(ratio, round2) * 100, '%', sep = '')]
    return(mdt[])
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ',',
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data <- read_file(infile1)
    res <- countSigSnpPerGenderPheno(data, pSigLevel, 4)
    write_file(res, outfile1)
}

run()


