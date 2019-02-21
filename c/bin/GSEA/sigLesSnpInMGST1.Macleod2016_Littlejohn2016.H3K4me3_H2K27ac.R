# /usr/bin/rscript

## Description:
## This module checks if the significant GSEA core SNPs are prioritised in
## Macleod 2016 and Littlejohn 2016 around the MGST1 gene.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 03 February 2017
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <infile2> <outfile1> <outfile2> <outfile3> <outfile4>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]    # lesESnp
infile2  <- args[2]    # Macleod 2016
infile3  <- args[3]    # Littlejohn 2016
outfile1 <- args[4]    # 


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
geneName <- f22(outfile1, '[.]', 7)[1]
pSiglvl  <- as.numeric(f22(outfile1, '[.]', 6)[1])
pheno    <- f22(outfile1, '[.]', 5)[1]
db       <- f22(outfile1, '[.]', 4)[1]


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

rawProcess1 <- function(data1, Pheno, Db, pSiglvl) {
    data1[, `:=` (pheno = Pheno, db = Db)]
    dat <- data1[pmeta <= pSiglvl]
    return(dat)
}

rawProcess2 <- function(data2, geneName) {
    dat <- data2[Gene_ID == geneName]
    snpNames <- paste('Chr', unlist(strsplit(dat$VariantPsition, ';')), sep = '')
    return(snpNames)
}

rawProcess3 <- function(data3) {
    res <- grep('^Chr', gsub('_', ':', data3$VariantName), value = TRUE)
    return(res)
}

getSep <- function(file) {
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
    Sep <- getSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- fread(infile2, header = TRUE)
    data3 <- fread(infile3, header = TRUE)

    data <- rawProcess1(data1, pheno, db, pSiglvl)
    hSnpNames <- unique(c(rawProcess2(data2, geneName), rawProcess3(data3)))   # highlighted SNP names 

    res <- data[lesESnpName %in% hSnpNames]
    write_file(res, outfile1)
}

run()



