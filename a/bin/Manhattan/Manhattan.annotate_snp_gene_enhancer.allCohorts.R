# /usr/bin/rscript

## Description:
## This script plots a manhattan plot over GWAS results across NA cohorts.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 18 July 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf> <gender> <pheno> <infile> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]  # GWAS result 
infile2  <- args[2]  # SNP annotation 
infile3  <- args[3]  # gene annotation 
infile4  <- args[4]  # putative enhancer sites 
outfile1 <- args[5]
# infile1  <- '/tmp//gender_pheno_chr_pos_effect_pvalue.tab'
# infile2  <- '/tmp//pheno_chr_pos_geneName.2016-01-19.a.tab'
# infile3  <- '/tmp//pheno_chr_start_end_geneName.2016-01-19.a.tab'
# infile4  <- '/tmp//pheno_chr_start_end_source.2016-01-19.a.tab'
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
logBase      <- as.numeric(f22(outfile1, '[.]', 6)[1]) 
nc           <- as.numeric(f22(outfile1, '[.]', 5)[1])
subPlotIndex <- f22(outfile1, '[.]', 4)[1]
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

getManhattanPos <- function(data1, Gender, Pheno) {
    data <- data1[gender == Gender & pheno == Pheno]
    n <- nchar(as.character(max(data$pos)))
    data[, manhattanPos := paste(sprintf("%02d", chr), sprintf(paste("%", n, "d", sep = ''), pos), sep = ':')]
    return(data)
}

log_transform <- function(pvalue, logBase) {
    return(-log(pvalue, logBase))
}

plot_Manhattan_allCohorts <- function(data1, Gender, Pheno, base, verysignifline, genomewideline, suggestiveline) {
    data <- getManhattanPos(data1, Gender, Pheno)
    Title <- paste(Gender, Pheno, sep =', ')
    Xlab  <- 'chromosome'
    Ylab  <- '-log10(p-value)'
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = manhattanPos, y = log_transform(pvalue, base), colour = as.factor(chr)), size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(verysignifline, base), colour = as.factor(verysignifline)), linetype = 4, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(genomewideline, base)), colour = as.factor(genomewideline), linetype = 2, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(suggestiveline, base)), colour = as.factor(suggestiveline), linetype = 1, size = 0.5) +
         scale_x_discrete(breaks = as.numeric(levels(as.factor(data$chr))), labels = as.numeric(levels(as.factor(data$chr)))) +
         scale_colour_discrete(guide = "none") +
         scale_alpha(guide = "none") +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 35),
               axis.title.x = element_text(size = 30),
               axis.title.y = element_text(size = 30),
               legend.text  = element_text(size = 25),
               legend.title = element_text(size = 25),
               axis.text.x  = element_text(size = 25),
               axis.text.y  = element_text(size = 25))
    return(g)
}

getLoopOverList <- function(data1) {
    genders <- unique(data1$gender)
    phenos  <- unique(data1$pheno)
    dt <- setDT(expand.grid(genders, phenos))
    setnames(dt, c('gender', 'pheno'))
    return(dt[])
}

draw_plots <- function(data1, loopList, subPlotIndex, nc, logBase, verysignifline, genomewideline, suggestiveline, outfile) {
    Title    <- 'Manhattan plots'
    Subtitle <- paste(paste('verysignifline', format(verysignifline, scientific = TRUE), sep = ':'),
                      paste('genomewideline', format(genomewideline, scientific = TRUE), sep = ':'),
                      paste('suggestiveline', format(suggestiveline, scientific = TRUE), sep = ':'), sep = '; ')
    nr <- nrow(loopList)/(nc-1) + 2
    Width  <- 4800 * nc
    Height <- 960 * (nr-2) + 960/10*3
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc + 1, heights = unit(c(2, 1, rep(10, nr-2)), rep("inches", nr)), widths = unit(c(1, 10), rep("inches", nc)))))
    grid.text(Title, gp = gpar(fontsize = 60, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:nc))
    grid.text(Subtitle, gp = gpar(fontsize = 50, fontface = "bold"), vp = viewport(layout.pos.row = 2, layout.pos.col = 1:nc))
    grid.text(subPlotIndex, gp = gpar(fontsize = 50, fontface = "bold"), vp = viewport(layout.pos.row = 1:nr, layout.pos.col = 1))
    for (i in seq(3, nr)) {
        index <- i - 2
        g <- plot_Manhattan_allCohorts(data1, loopList[index]$gender, loopList[index]$pheno, logBase, verysignifline, genomewideline, suggestiveline)
        print(g, vp = viewport(layout.pos.row = i, layout.pos.col = nc))
    }
    dev.off()
}

run <- function() {
    data1 <- read_file(infile1)
    loopList <- getLoopOverList(data1)
    draw_plots(data1, loopList, subPlotIndex, nc, logBase, verysignifline, genomewideline, suggestiveline, outfile1)
}

run()


