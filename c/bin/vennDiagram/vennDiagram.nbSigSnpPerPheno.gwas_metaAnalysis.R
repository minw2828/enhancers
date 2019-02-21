# /usr/bin/rscript

## Description:
## This script plots a venn diagram to visualise the overlap of significant SNPs 
## between GWAS and meta-analysis.
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
infile2  <- args[2]  # meta-analysis result 
outfile1 <- args[3]
# infile1 <- '/tmp//gender_pheno_chrN_snpName_effect_pvalue.2016-01-19.c.csv'
# infile2 <- '/tmp//chrN_snpName_pheno_snpeffectmeta_pmeta.2016-01-19.c.csv'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.nbSigSnpPerGenderPheno/gender_pheno_nbSigSnp.10e-08.2016-01-19.c.csv'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pSigLevel <- f22(outfile1, '[.]', 4)[1]
analyses  <- c('GWAS', 'Meta-analysis')


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getSigSnpIdGwas <- function(data1, pSigLevel) {
    pSigLevel <- as.numeric(pSigLevel)
    dat1 <- data1[pvalue <= pSigLevel]
    dat1[, sigSnpId := paste(gender, pheno, snpName, sep = ':')]
    return(dat1)
}

getSigSnpIdMeta <- function(data2, pSigLevel) {
    pSigLevel <- as.numeric(pSigLevel)
    dat2.sig  <- data2[pmeta <= pSigLevel]
    dat2.sig2 <- copy(dat2.sig)
    dat2.sig[, gender := 'bull']
    dat2.sig2[, gender := 'cow']
    dat2 <- rbindlist(list(dat2.sig, dat2.sig2))
    dat2[, sigSnpId := paste(gender, pheno, snpName, sep = ':')]
    return(dat2)
}

venn_diagram <- function(data1.sig, data2.sig, pSigLevel, analyses) {
    Title    <- paste('All: GWAS and Meta-analysis significant variants')
    Subtitle <- paste('significance level', pSigLevel, sep = ': ')
    venn <- venn.diagram(x = list(GWAS = data1.sig$sigSnpId, META = data2.sig$sigSnpId), 
                         filename = NULL, height = 3000, width = 3000, resolution = 500,
                         imagetype = "png", units = "px", compression = "lzw", na = "stop", 
                         main = Title, sub = Subtitle, 
                         main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif", 
#                         main.col = "black", main.cex = 3, main.just = c(0.5, 1), 
                         main.col = "black", main.cex = 4.8, main.just = c(0.5, 1),
                         sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif", 
#                         sub.col = "black", sub.cex = 2, sub.just = c(0.5, 1), 
                         sub.col = "black", sub.cex = 4.8, sub.just = c(0.5, 1),
                         category.names = analyses, cat.col = c("royalblue1", "hotpink1"), 
#                         cat.pos = c(-50, 50), cat.dist = rep(0.025, 2), cat.cex = rep(2.5, 2),
                         cat.pos = c(-40, 40), cat.dist = rep(0.025, 2), cat.cex = rep(5, 2),
                         force.unique = TRUE, print.mode = "raw", sigdigs = 3, 
                         direct.area = FALSE, area.vector = 0, hyper.test = FALSE, 
                         total.population = NULL, lower.tail = TRUE,
                         col = c("skyblue1", "pink1"), fill = c("skyblue1", "pink1"), alpha = rep(0.4, 2), 
#                         cex = 2.5, scaled = TRUE)
                         cex = 6, scaled = TRUE)
    return(venn)
}

draw_plots <- function(g, outfile) {
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    nr <- 1;     nc <- 1
    rhText <- 1; rhPlot <- 5
    nText <- 0;  nPlot <- 1
    Width  <- 1200 * nc
    Height <- 960 * (rhText * rhText + rhPlot * nPlot)/rhPlot
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    plot.new() # start new page
    gl <- grid.layout(nrow = nr, ncol = nc, heights = unit(c(rep(rhText, nText), rep(rhPlot, nPlot)), "null")) # setup layout
    # grid.show.layout(gl)
    vp.1 <- viewport(layout.pos.row = 1, layout.pos.col = 1) # setup viewports
    pushViewport(viewport(layout = gl)) # init layout
    pushViewport(vp.1) # access the first position
    grid.draw(g)
    popViewport() # done with the second viewport
    dev.off()
}


run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    data1.sig <- getSigSnpIdGwas(data1, pSigLevel)
    data2.sig <- getSigSnpIdMeta(data2, pSigLevel) 
    g <- venn_diagram(data1.sig, data2.sig, pSigLevel, analyses)
    draw_plots(g, outfile1)
}

run()


