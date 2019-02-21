# /usr/bin/rscript

## Description:
## This module selects genes of interest.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 19 May 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <infile2> <outfile1> <outfile2> <outfile3>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile  <- args[3]
# infile1  <- '/group/dairy/Min/geno2pheno/data/qtl/ncbi.cattleQTLmilk.20160517.tab'
# infile2  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/enriched/GSEAenriched.gender_pheno_chrN_db.0.0.01.10000.2016-01-19.c.csv'
# outfile  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/genesFound/ncbi.cattleQTLmilk.20160517.1000.tab'

# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(infile2, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(infile2, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
outpath  <- file.path(path_pre, 'analyses', date, psf, 'out')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

# Get global variables #
distance <- as.numeric(f22(outfile, '[.]', 2)[1])
infile3  <- 'output.geneName_geneChr_geneStart_geneEnd_sigSNPs'

# analysis #
rawProcess <- function(data) {
    return(data[!is.na(start_position_on_the_genomic_accession)])
}

get_targets <- function(data1, data2) {
    targets <- data2[Chromosome %in% data1$chromosome]
    targets[, files := file.path(outpath, 'GSEA', Gender, Phenotype, sprintf("%02d", Chromosome), Database, paste(infile3, Gender, Phenotype, sprintf("%02d", Chromosome), Database, distance, date, psf, 'tab', sep = '.'))]
    return(targets)
}

check_geneFound <- function(data1, target) {
    geneNames <- data1[chromosome == target$Chromosome]$Symbol
    data <- setnames(fread(target$file), c('geneName', 'chromosome', 'geneStart', 'geneEnd', 'snpNearby'))
    res <- sapply(geneNames, function(x) length(grep(x, data$geneName)))
    result <- data.table(Symbol = geneNames,
                         Gender = target$Gender,
                         Phenotype = target$Phenotype,
                         Chromosome = target$Chromosome,
                         Database = target$Database,
                         boolFound = res) # 0: not found; !0: found
    return(result)
}

result <- function(data1, data2) {
    setkey(data1, Symbol)
    setkey(data2, Symbol)
    return(merge(data1, data2)[, .(Symbol, description, start_position_on_the_genomic_accession, end_position_on_the_genomic_accession, Gender, Phenotype, Chromosome, Database, boolFound)])
}

## write file ##
write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
}

## run ##
run <- function() {
    data1 <- rawProcess(fread(infile1))
    data2 <- fread(infile2)
    targets <- get_targets(data1, data2)
    res <- rbindlist(lapply(seq(nrow(targets)), function(x) check_geneFound(data1, targets[x])))
    content <- result(data1, res)
    write_file(content, outfile)
}

run()



