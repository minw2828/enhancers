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
pTtest     <- as.numeric(f22(outfile1, '[.]', 5)[1])
ORefseq    <- f22(outfile1, '[.]', 4)[1]


# get supporting functions #
#Task <- f21(f22(program, '/', 1), '[.]', 1)


# set seed #
set.seed(1234)    # Need to set seed when using caret package


# analysis #
getResult <- function(data, pTtest) {
    data[, boo := ifelse(pvalue <= pTtest, TRUE, FALSE)]
    return(data[])
}

run <- function() {
    data1 <- read_file(infile1)

    data <- getResult(data1, pTtest)
    write_file(data, outfile1, FALSE)
}

run()



