# /usr/bin/rscript

## Description:
## This module calculates the variance of SNP effect.
##
## Reference:
## 2014, Fragomeni B. de O. et al., Changes in variance explained by top SNP 
## windows over generations for three traits in broiler chicken, Front Genet
## 2010, Zhang Z. et al., Best Linear Unbiased Prediction of Genomic Breeding 
## Values Using a Trait-Specific Marker-Derived Relationship Matrix, PLoS One
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 16 March 2015
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <infile2> <outfile> <outfile2> <outfile3>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/bull/FY/01/chrom01.ps'
# infile2  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/bull/FY/01/VISTA/snpName_freqAllele0_freqAllele1_freqAlleleMinor_freqGeno0_freqGeno1_freqGeno2.bull.FY.01.VISTA.2016-01-19.c.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/bull/FY/01/VISTA/snpName_effect_pvalue_variance.bull.FY.01.VISTA.2016-01-19.c.tab'

# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

# Get global variables #
gender <- f21(outfile1, '[.]', 2)
pheno  <- f21(outfile1, '[.]', 3)
chrN   <- as.numeric(f21(outfile1, '[.]', 4))
db     <- f21(outfile1, '[.]', 5)

## read file and give it colname ##
read_file <- function(infile, cns) {
    return(setnames(fread(infile, header = FALSE), cns))
}

calculate_varianceOfSnpEffect <- function(data1, data2) {
    setkey(data1, snpName)
    setkey(data2, snpName)
    data <- merge(data1, data2)
    data[, variance := effect ** 2 * 2 * freqAllele0 * freqAllele1]
    return(data[,.(snpName, effect, pvalue, variance)])
}

## write file ##
write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
}

# speed check 
# system.time(.fun)

## run ##
run <- function() {
    data1 <- read_file(infile1, c('snpName', 'effect', 'pvalue'))
    data2 <- read_file(infile2, unlist(strsplit(f22(f21(infile2, '[.]', 1), '/', 1), '_')))
    write_file(calculate_varianceOfSnpEffect(data1, data2), outfile1)
}

run()

