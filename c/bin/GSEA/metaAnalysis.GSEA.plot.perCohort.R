# /usr/bin/rscript

## Description:
## This module plots up results from GSEA based on meta-analysis result.
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
## Rscript <module_name> <infile1> <infile2> <infile3> <infile4> <outfile1> <outfile2> <outfile3> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]  # GSEA cal result 
infile2  <- args[2]  # ES
infile3  <- args[3]  # esNulls
infile4  <- args[4]  # EsRank
outfile1 <- args[5]
outfile2 <- args[6]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/FY/VISTA/pheno_chrN_snpName_snpeffectmeta_varmeta_zmeta_pmeta_minusLog10Pmeta_boolenhancer_pHitMiss_accSum.FY.VISTA.2016-01-19.c.tab'
# infile2  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/FY/VISTA/ES.FY.VISTA.10000.2016-01-19.c.tab'
# infile3  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/FY/VISTA/esNull.FY.VISTA.10000.2016-01-19.c.tab'
# infile4  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/FY/VISTA/EsRank.FY.VISTA.10000.2016-01-19.c.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/FY/VISTA/plot12.FY.VISTA.10000.2016-01-19.c.png'
# outfile2 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/FY/VISTA/plot123456.FY.VISTA.10000.2016-01-19.c.png'


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pheno  <- f22(outfile1, '[.]', 6)[1]
db     <- f22(outfile1, '[.]', 5)[1]
ntimes <- as.integer(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file1 <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

read_file2 <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- as.data.table(read.table(infile))
    setnames(data, cns)
    return(data)
}

get_Theme <- function() {
    # element_text(family, face, colour, size)
    return(theme(plot.title = element_text(face = "bold", size = 25),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20),
                 axis.text.x= element_text(size = 20),
                 axis.text.y = element_text(size = 20),
                 axis.text.x= element_text(size = 20),
                 axis.text.y = element_text(size = 20)))
}

plot_running_sum <- function(data) {
    Title <- 'Running Sum: sorted -log10(pmeta) (largest -> smallest)'
    Xlab  <- "-log10(pmeta) Rank"
    Ylab  <- "Running Enrichment Score"
    legend_name_fill <- "boolean enhancer"
    data[, rIndex := data[,.I]]
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = rIndex, y = accSum), colour = "#009900", size = 1) +
         geom_hline(mapping = aes(yintercept = 0), colour = "#009900", linetype = 1, size = 1) +
         geom_rect(mapping = aes(xmin = rIndex, xmax = rIndex, ymin = max(accSum)+max(accSum)/10, ymax = Inf, fill = as.factor(boolenhancer)), size = 1) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_histogram_esNull_ES <- function(ESs, ES, EsRank) {
    Title <- 'Histogram: Enrichment Scores ES(NULL).  X-intercept: ES(S)'
    bin_size <- 0.01
    Xlab  <- paste(paste("ES(NULL);   bin size:", bin_size), ';   ', paste(ntimes, 'permutations'), sep = "")
    Ylab  <- "Frequency"
    legend_name_fill     <- paste("ES(NULL)", Ylab)
    legend_name_linetype <- "ES(S)"
    g <- ggplot(data = ESs) +
         geom_histogram(mapping = aes(esNull, fill = ..count..), binwidth = bin_size) +
         scale_fill_gradient(name = legend_name_fill, low = "#99FF33", high = "#009900") +
         geom_vline(mapping = aes(xintercept = ES, linetype = 'k'), colour = "#003300", size = 1) +
         geom_text(mapping = aes(x = ES, y = 0, label = paste('ES(S). Rank', EsRank, sep =': ')), alpha = 0.4, size = 8, colour = "#003300", hjust = -1, vjust = 1.5, angle = 90) +
         scale_linetype_manual(legend_name_linetype, values = c('k' = 2), labels = c(round(ES, 2))) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_histogram_minusLog10Pmeta_boolenhancer <- function(data) {
    Title <- 'Histogram: Distribution of -log10(pmeta)'
    bin_size <- 10
    Xlab  <- paste("-log10(pmeta);  bin size:", bin_size)
    Ylab  <- "Frequency"
    legend_name_fill <- "boolean enhancer"
    g <- ggplot(data = data) +
         geom_histogram(mapping = aes(minusLog10Pmeta, fill = boolenhancer), binwidth = bin_size) +
         facet_wrap(~ boolenhancer, ncol = 1) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_boxplot_minusLog10PmetaRank_boolenhancer <- function(data) {
    Title <- 'Boxplot: -log10(pmeta) Rank'
    Xlab  <- "boolean enhancer"
    Ylab  <- "-log10(pmeta) Rank"
    legend_name_fill <- "boolean enhancer"
    data[, rIndex := data[,.I]]
    g <- ggplot(data = data) +
         geom_boxplot(mapping = aes(boolenhancer, rIndex, fill = boolenhancer)) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_snpeffectmeta_minusLog10Pmeta <- function(data) {
    Title <- 'Scatterplot: effectmeta v.s. -log10(pmeta)'
    Xlab  <- "effectmeta"
    Ylab  <- "-log10(pmeta)"
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = snpeffectmeta, y = minusLog10Pmeta), colour = "#009900", size = 1) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_minusLog10PmetaRank_pHitMiss <- function(data) {
    Title <- 'Scatterplot: pHits/pMiss v.s. -log10(pmeta) Rank'
    Xlab  <- "-log10(pmeta) Rank"
    Ylab  <- "pHitMiss"
    data[, rIndex := data[,.I]]
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = rIndex, y = pHitMiss, colour = boolenhancer), size = 1) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

draw_plots12 <- function(g1, g2, outfile) {
    nr <- 3
    nc <- 1
    Height_unit <- c(1, 5, 5)
    Height <- 960/5*sum(Height_unit)
    Width  <- 2400*nc
    Text <- paste(pheno, db, sep = ', ')
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(Height_unit, "null"))))
    grid.text(Text, gp = gpar(fontsize = 40), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(g1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(g2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    dev.off()
}

draw_plots123456 <- function(g1, g2, g3, g4, g5, g6, outfile) {
    nr <- 5
    nc <- 5
    Height_unit <- c(1, 5, 5, 5, 5)
    Height <- 960/5*sum(Height_unit)
    Width  <- 1200*nc
    Text <- paste(pheno, db, sep = ', ')
    # plot distbution ratio here is: Width : Height = (2400/2*3) : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(Height_unit, "null"))))
    grid.text(Text, gp = gpar(fontsize = 100), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
    print(g1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:4))
    print(g2, vp = viewport(layout.pos.row = 2, layout.pos.col = 5))
    print(g3, vp = viewport(layout.pos.row = 3, layout.pos.col = 5))
    print(g4, vp = viewport(layout.pos.row = 4, layout.pos.col = 5))
    print(g5, vp = viewport(layout.pos.row = 5, layout.pos.col = 5))
    print(g6, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:4))
    dev.off()
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
}

run <- function() {
    data1 <- read_file1(infile1)  # GSEA cal result
    data2 <- read_file1(infile2)  # ES 
    data3 <- read_file1(infile3)  # esNulls
    data4 <- read_file2(infile4)  # EsRank
    g1 <- plot_running_sum(data1)
    g2 <- plot_histogram_esNull_ES(data3, data2$ES, data4$EsRank)
    g3 <- plot_histogram_minusLog10Pmeta_boolenhancer(data1)
    g4 <- plot_boxplot_minusLog10PmetaRank_boolenhancer(data1)
    g5 <- plot_snpeffectmeta_minusLog10Pmeta(data1)
    g6 <- plot_minusLog10PmetaRank_pHitMiss(data1)
    draw_plots12(g1, g2, outfile1)
    draw_plots123456(g1, g2, g3, g4, g5, g6, outfile2)
}

run()



