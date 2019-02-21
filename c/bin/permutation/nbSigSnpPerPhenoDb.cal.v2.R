# /usr/bin/rscript

## Description:
## On 14 September 2017, MG, AC, TPH and MW had a skype call to BH. After the
## meeting, TPH requests MW to perform the permutation test using MG's sliding
## method, in replacement of the random selection permutation method that was
## used in the paper.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 15 September 2017
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]    # TAD info
infile2  <- args[2]    # metaAnalysis 
infile3  <- args[3]    # db SNPs
outfile1 <- args[4]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
source(file.path(binpath, 'getGenome.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables # 
pheno   <- f22(outfile1, '/', 4)[1]
db      <- f22(outfile1, '/', 3)[1]
Chr     <- as.numeric(f22(outfile1, '/', 2)[1])
pSig    <- as.numeric(f22(outfile1, '[.]', 6)[1])
ntimes  <- as.numeric(f22(outfile1, '[.]', 5)[1])
ORefseq <- f22(outfile1, '[.]', 4)[1]


# set seeds
set.seed(1234)


# analysis #
getCoverage <- function(gr2, gr3, pSig) {
    sgr2 <- reduce(gr2[gr2$pmeta <= pSig])
    cr <- GenomicRanges :: coverage(GRangesList(sgr2, gr3))
    scr <- slice(cr, lower = 2)
    res <- sum(as.numeric(unlist(width(scr))))
    return(res)
}

getRoundBack <- function(gr, moveN, genomeLen, Chr) {
    ngr <- GenomicRanges::shift(gr, moveN)
    tmp1 <- ngr[end(ngr) <= genomeLen]
    tmp2 <- ngr[start(ngr) <= genomeLen & end(ngr) > genomeLen]
    tmp3 <- ngr[start(ngr) > genomeLen]
    if (length(tmp2) != 0) {
        tmp2A <- copy(tmp2)
        end(tmp2A) <- genomeLen
        tmp2B <- GRanges(seqnames = Chr, IRanges(start = 1, end = end(tmp2)-genomeLen))
    } else {
        tmp2A <- tmp2B <- GRanges()
    }
    if (length(tmp3) != 0) {
        start(tmp3) <- start(tmp3) - genomeLen
        end(tmp3) <- end(tmp3) - genomeLen
    } else {
        tmp3 <- GRanges()
    }
    res <- c(tmp1, tmp2A, tmp2B, tmp3)
    return(res)
}

shiftTadbByChr <- function(gr, genomeLength) {
    moveN <- sample(1:genomeLength, 1)
    res <- getRoundBack(gr, moveN, genomeLength, Chr)
    return(res)
}

getNullCoverage <- function(gr1, gr2, genomeLength, pSig) {
    cgr1 <- copy(gr1)
    isCircular(cgr1) <- rep(TRUE, length(isCircular(cgr1)))
    ns <- as.list(levels(seqnames(cgr1)))
    tmp <- endoapply(ns, function(x) shiftTadbByChr(cgr1, genomeLength))
    ngr1 <- do.call(c, tmp)
    isCircular(ngr1) <- rep(NA, length(isCircular(ngr1)))
    isCircular(gr2) <- rep(NA, length(isCircular(gr2)))
    res <- getCoverage(ngr1, gr2, pSig)
    return(res)
}

calculate_permutations <- function(gr1, gr2, genomeLength, pSig, n) {
    res <- replicate(n, getNullCoverage(gr1, gr2, genomeLength, pSig), simplify = "vector")
    return(res)
}

getOriPermu <- function(gr2, gr3, genomeLength, pSig, ntimes) {
    ori <- getCoverage(gr2, gr3, pSig)
    permu <- calculate_permutations(gr2, gr3, genomeLength, pSig, ntimes)
    return(c(ori, permu))
}

run <- function() {
    data1 <- read_file(infile1) 
    data2 <- read_file(infile2) 
    data3 <- read_file(infile3)

    ah <- AnnotationHub()
    genome <- getGenome(data1, NULL, ORefseq)
    genomeLength <- as.numeric(seqlengths(genome)[Chr])

    gr2 <- getGranges(data2) # must be this way, do reduce later
    gr3 <- reduce(getGranges(data3))

    res <- getOriPermu(gr2, gr3, genomeLength, pSig, ntimes)
    write_file(res, outfile1, FALSE)
}

run()



