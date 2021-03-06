# /usr/bin/rscript

## Description:
## On 14 September 2017, MG, AC, TPH and MW had a skype call to BH. MG raises
## a concern on the permutation test that MW performed on her first paper.
## MG's concern was that whether the enrichment signal from H3K4me3 regions
## was due to LD to mutation in gene, or true enrichment within H3K4me3 regions.
## To address the concern, MG proposes the following analyses that will be
## performed in this folder.
##
## This script makes a table according to MG's suggestions with the following
## columns:
## - Column 1: Gene name
## - Column 2: ratio of the upstream that are covered by this histone mark
## - Column 3: number of QTLs/eQTLs falling inside that histone mark
## - Column 4: number of QTLS/eQTLs falling outside that histone mark
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 14 September 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <ofnpattern>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1   <- args[index+1]    # 
outfile1  <- args[index+2]    # 


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
outpath  <- file.path(path_pre, 'analyses', date, psf, 'out')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
source(file.path(binpath, 'getData1234.R'))
source(file.path(binpath, 'getGTF.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables #
cov1    <- as.numeric(f22(outfile1, '[.]', 7)[1])
cov2    <- as.numeric(f22(outfile1, '[.]', 6)[1])
pSig    <- as.numeric(f22(outfile1, '[.]', 5)[1])
ORefseq <- f22(outfile1, '[.]', 4)[1]


# get supporting functions #
#Task <- f21(f22(program, '/', 1), '[.]', 1)


# set seed #
set.seed(1234)    # Need to set seed when using caret package


# analysis #
rawProcess1 <- function(data1, cov1, cov2) {
    data1[, boo := ifelse(nbQltInDb >  nbQltOutDb, 1, 
                   ifelse(nbQltInDb == nbQltOutDb, 0, -1))]
    dat1 <- data1[ratioOfProInDb >= cov1 & ratioOfProInDb <= cov2]
    return(dat1)
}

run <- function() {
    data1 <- read_file(infile1)
    data1[, boo := ifelse(nbQltInDb >  nbQltOutDb, 1,
                   ifelse(nbQltInDb == nbQltOutDb, 0, -1))]
    dt <- setDT(ddply(data1, .(pheno, db), summarize, paste(paste(names(table(boo)), table(boo), sep = ':'), collapse = '; ')))
 
    dat1 <- rawProcess1(data1, cov1, cov2)
    nrow(dat1)
    table(dat1$boo)

    res <- getResult(data2, data3, gtf, pSig, pDis1, pDis2)
    write_file(res, outfile1, FALSE)
}

run()



