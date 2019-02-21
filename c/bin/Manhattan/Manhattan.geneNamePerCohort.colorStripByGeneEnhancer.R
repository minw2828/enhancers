# /usr/bin/rscript

## Description:
## This script plots a manhattan plot over meta-analysis results.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 20 July 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]  # gene info 
infile2  <- args[2]  # meta-analysis result on the same chromosome as gene 
infile3  <- args[3]  # enhancer regions on the same chromosome as gene 
outfile1 <- args[4]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
geneName     <- f22(outfile1, '[.]', 8)[1]
logBase      <- as.numeric(f22(outfile1, '[.]', 7)[1]) 
Width        <- as.numeric(f22(outfile1, '[.]', 6)[1])
Height       <- as.numeric(f22(outfile1, '[.]', 5)[1])
nc           <- as.numeric(f22(outfile1, '[.]', 4)[1])
verysignifline <- 10e-20
genomewideline <- 10e-8
suggestiveline <- 10e-5


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

rawProcess2 <- function(data2, logBase) {
    rpmeta <- head(unique(data2[order(pmeta)]$pmeta), 2)[2]
    data2[pmeta == 0, pmeta := rpmeta]
    data2[, logP := -log(pmeta, logBase)]
}

plot_Manhattan <- function(data1, data2, data3, verysignifline, genomewideline, suggestiveline) {
    Title <- paste(unique(data2$pheno), data1$geneName, sep = ': ')
    Xlab  <- paste('\nchromosome', paste(unique(data1$chr), '\n', sep = ''), sep = ': ')
    Ylab  <- '-log10(p-value)\n'
    legend_name_colour <- "SNP P-value & significance threshold: "
    legend_name_fill <- "gene & enhancer: "
    g <- ggplot() +
         geom_point(mapping = aes(x = pos, y = logP), data = data2, size = 15) +
         geom_hline(mapping = aes(yintercept = verysignifline, colour = as.factor(verysignifline)), linetype = 4, size = 10, alpha = 0.8, show.legend = TRUE) +
         geom_hline(mapping = aes(yintercept = genomewideline, colour = as.factor(genomewideline)), linetype = 2, size = 10, alpha = 0.8, show.legend = TRUE) +
         geom_hline(mapping = aes(yintercept = suggestiveline, colour = as.factor(suggestiveline)), linetype = 1, size = 10, alpha = 0.8, show.legend = TRUE) +
         geom_rect(mapping = aes(xmin = start, ymin = -5, xmax = end, ymax = Inf, fill = geneName), data = data1, alpha = 0.4, show.legend = TRUE) +
         geom_rect(mapping = aes(xmin = start, ymin = -3, xmax = end, ymax = Inf, fill = db), data = data3, alpha = 0.4, show.legend = TRUE) +
         scale_alpha(guide = "none") +
         labs(title = Title, x = Xlab, y = Ylab) +
         guides(colour = guide_legend(title = legend_name_colour, nrow = 1, keywidth = 3, keyheight = 3, default.unit = "inch", override.aes = list(fill = NA)),
                fill   = guide_legend(title = legend_name_fill, nrow = 1, override.aes = list(colour = NA))) +
         theme(plot.title   = element_text(face = "bold", size = 200),
               axis.title.x = element_text(size = 180), axis.title.y = element_text(size = 180),
               axis.text.x  = element_text(size = 180), axis.text.y  = element_text(size = 180),
               legend.title = element_text(size = 180), legend.text  = element_text(size = 180), 
               legend.key.size = unit(5.5, 'cm'),    legend.position = "bottom",
               plot.margin = unit(c(5,5,5,5), "cm"))
    return(g)
}

run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    data3 <- read_file(infile3)

    rawProcess2(data2, logBase)

    geneNames <- unique(data1$geneName)
    dis <- 10000
    dat2 <- rbindlist(lapply(geneNames, function(x) getSnpAroundGene(data1, data2, x, dis)))
    dat3 <- rbindlist(lapply(geneNames, function(x) getEnhancerAroundGene(data1, data3, x, dis)))

    g <- plot_Manhattan(data1, dat2, dat3, -log(verysignifline, logBase), -log(genomewideline, logBase), -log(suggestiveline, logBase))
    f131(outfile1, g, Width, Height)
}

run()


