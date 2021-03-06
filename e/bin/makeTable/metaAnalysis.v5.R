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
## All the genes are output from the results, regardless of whether there is any
##   GWAS hit inside the gene or not.
## The QTL/eQTL data that is being tested for enrichment is listed as function name.
##
## Version 5 is developed from Version 4:
## - Get genomic regions that have no annotation at all, to an extend that sequences
##   on both strands do not have any annotations at all. See if histone within those
##   regions are getting more variants or not.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 28 September 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <ofnpattern>


# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1   <- args[index+1]    # TAD info
infile2   <- args[index+2]    # metaAnalysis 
infile3   <- args[index+3]    # database 
infile4   <- args[index+4]    # NCBI annotation
outfile1  <- args[index+5]    # 


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
source(file.path(binpath, 'getGenome.R'))
source(file.path(binpath, 'getGap.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables #
tmp1 <- as.numeric(f22(outfile1, '/', 5)[1])
tmp2 <- as.numeric(f22(outfile1, '/', 4)[1])
pDis1     <- ifelse(is.na(tmp1), NA, tmp1)
pDis2     <- ifelse(is.na(tmp2), NA, tmp2)
pheno    <- f22(outfile1, '/', 3)[1]
db       <- f22(outfile1, '/', 2)[1]
pSigCTCF <- ifelse(db == 'CTCF', as.numeric(f22(outfile1, '[.]', 7)[1]), NA)
ms       <- ifelse(db == 'CTCF', as.numeric(f22(outfile1, '[.]', 6)[1]), NA)
pSigMeta <- as.numeric(f22(outfile1, '[.]', 5)[1])
ORefseq  <- f22(outfile1, '[.]', 4)[1]


# get supporting functions #
Task <- f22(program, '/', 2)[1]
source(file.path(binpath, Task, 'getDbInPromoter.R'))


# set seed #
set.seed(1234)    # Need to set seed when using caret package


# analysis #
getRgr34 <- function(data3, data4, db, pSigCTCF, ms) {
    if (db == 'CTCF') {
        res <- reduce(getGranges(getData34(data3, data4, pSigCTCF, ms)))
    } else {
        res <- reduce(getGranges(data3))
    }
    return(res)
}

getSliced <- function(pro, chr, sgr3, nCov) {
    if (length(sgr3) != 0) {
        cr <- coverage(GRangesList(reduce(pro), reduce(sgr3)))
        scr <- cr[[chr]]
        sscr <- slice(scr, lower = nCov)
        if (length(sscr) != 0)  {
            res <- GRanges(seqnames = chr, ranges(sscr))
        } else {
            res <- GRanges()
        }
    } else {
        res <- GRanges()
    }
    return(res)
}

getDbInPro <- function(pro, rgr3, db, chr, nCov) {
    sgr3 <- rgr3[seqnames(rgr3) == chr]
    if (db == 'CTCF') {
        sgr31 <- sgr3[strand(sgr3) == '+']
        sgr32 <- sgr3[strand(sgr3) == '-']
        tmp1 <- getSliced(pro, chr, sgr31, nCov)
        tmp2 <- getSliced(pro, chr, sgr32, nCov)
        res <- c(tmp1, tmp2)
    } else {
        res <- getSliced(pro, chr, sgr3, nCov)
    } 
    return(res)
}

get1Record <- function(sgr2, rgr3, gr, genome, pDis1, pDis2, db, nCov, x) {
    print(x)
    gene <- gr[x]
    chr <- as.character(seqnames(gene))
    genomeLen <- as.numeric(seqlengths(genome)[chr])
    his <- getDbInPro(gene, rgr3, db, chr, nCov)
    tmp1 <- data.table(intergenicName = paste(seqnames(gr[x]), start(gr[x]), end(gr[x]), strand(gr[x]), sep = '-'),
                       ratioOfProInDb = getRatioOfProInDb(gene, his))
    tmp2 <- getNbQtlInOutDb(sgr2, gene, his)
    res <- cbind(tmp1, tmp2)
    return(res)
}

getResult <- function(data2, rgr3, gr, genome, pSig, pDis1, pDis2, db, nCov) {
    dat2 <- data2[pmeta <= pSig]
    sgr2 <- getGranges(dat2)
    ns <- seq(length(gr))
    res <- rbindlist(lapply(ns, function(x) get1Record(sgr2, rgr3, gr, genome, pDis1, pDis2, db, nCov, x)))
    return(res)
}

getGaps <- function(genome, Refseq, gr) {
    anno <- getSeqlens(genome, Refseq)
    seqlengths(gr) <- anno[names(anno) %in% names(seqlengths(gr))]
    tmp  <- gaps(gr)
    res  <- tmp[strand(tmp) == '*']
    return(res)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    data3 <- read_file(infile3)
    data4 <- fread(infile4)

    ah <- AnnotationHub()
    genome <- getGenome(data1, NULL, ORefseq) 
    gtf <- getGtfBySpeciesRefseq(data1, ah, ORefseq, c('gene', 'transcript', 'exon', 'CDS', 'start_codon', 'stop_codon', 'three_prime_utr', 'five_prime_utr'), c('gene_name', 'gene_id'))
    strand(gtf) <- '*'
    gr <- getGaps(genome, Refseq, gtf)   

    gapRatio <- sum(as.numeric(width(gr))) / sum(as.numeric(seqlengths(genome))) * 100
    answer <- paste('The strand inspecific gaps takes up', gapRatio, '% of the genome.')
    print(answer)

    rgr34 <- getRgr34(data3, data4, db, pSigCTCF, ms)

    res <- getResult(data2, rgr34, gr, genome, pSigMeta, pDis1, pDis2, db, 2)
    write_file(res, outfile1, FALSE)
}

run()



