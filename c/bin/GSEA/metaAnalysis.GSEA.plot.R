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
infile2  <- args[2]  # esNulls
infile3  <- args[3]  # EsRank
outfile1 <- args[4]
outfile2 <- args[5]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/snpName_pheno_chrN_snpeffectmeta_varmeta_zmeta_pmeta_minusLog10Pmeta_source_boolenhancer_pHitMiss_accSum.2016-01-19.c.tab'
# infile2  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/esNull_id_pheno_db.10000.2016-01-19.c.tab'
# infile3  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/EsRank_pheno_db.10000.2016-01-19.c.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/runningSum_esNullES.10000.2016-01-19.c.png'
# outfile2 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis.GSEA/runningSum_esNullES_variance_rankVariance_effectVariance_pHitMiss_lesRunningSum_lesPHitMiss.10000.2016-01-19.c.png'


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
ntimes <- as.integer(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getLoopOverList <- function(data) {
    genders <- unique(data$gender)
    phenos  <- unique(data$pheno)
    dbs     <- unique(data$source)
    dt <- setDT(expand.grid(genders, phenos, dbs))
    setnames(dt, c('gender', 'pheno', 'db'))
    dt[, pointColor := colorRampPalette(brewer.pal(8, "Set2"))(nrow(dt))]
    dt[, lineColor  := colorRampPalette(brewer.pal(8, "Set1"))(nrow(dt))]
    dt[, lowFill    := colorRampPalette(brewer.pal(8, "Dark2"))(nrow(dt))]
    dt[, highFill   := colorRampPalette(brewer.pal(8, "Dark2"))(nrow(dt))]
    return(dt[])
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

## generate running-sum statistic plot. Write to outfile ##
### ggplot2 - Easy way to mix multiple graphs on the same page - R software and data visualization
### http://www.sthda.com/english/wiki/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page-r-software-and-data-visualization
plot_running_sum <- function(data, pointColor) {
    Title <- paste(unique(data$pheno), unique(data$source))
#    Title <- 'Running Sum: sorted -log10(p-value) (largest -> smallest)'
    Xlab  <- "-log10(p-value) Rank"
    Ylab  <- "Running Enrichment Score"
    legend_name_fill <- "boolean enhancer"
    data[, rIndex := data[,.I]]
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = rIndex, y = accSum), colour = pointColor, size = 1) +
         geom_hline(mapping = aes(yintercept = 0), colour = pointColor, linetype = 1, size = 1) +
         geom_rect(mapping = aes(xmin = rIndex, xmax = rIndex, ymin = max(accSum)+max(accSum)/10, ymax = Inf, fill = as.factor(boolenhancer)), size = 1) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

## generate historgram of variance: (enhancer) & (enhancer & non-enhancer) ##
plot_histogram_variance_boolenhancer <- function(data) {
    Title <- 'Histogram: Distribution of -log10(p-value)'
    bin_size <- 100
    Xlab  <- paste("-log10(p-value);  bin size:", bin_size)
    Ylab  <- "Frequency"
    legend_name_fill <- "boolean enhancer"
    g <- ggplot(data = data) +
         geom_histogram(mapping = aes(variance, fill = boolenhancer), binwidth = bin_size) +
         facet_wrap(~ boolenhancer, ncol = 1) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

## generate historgram of variance: (enhancer) & (enhancer & non-enhancer) ##
plot_effect_variance <- function(data) {
    Title <- 'Scatterplot: effect v.s. -log10(p-value)'
    Xlab  <- "effect"
    Ylab  <- "-log10(p-value)"
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = effect, y = variance), colour = "#009900", size = 1) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

## plot histogram of enrichment scores esNull ##
plot_histogram_esNull_ES <- function(ESs, ES, ES.rank) {
    Title <- 'Histogram: Enrichment Scores ES(NULL).  X-intercept: ES(S)'
    bin_size <- 0.001
    Xlab  <- paste(paste("ES(NULL);   bin size:", bin_size), ';   ', paste(ntimes, 'permutations'), sep = "")
    Ylab  <- "Frequency"
    legend_name_fill     <- paste("ES(NULL)", Ylab)
    legend_name_linetype <- "ES(S)"
    g <- ggplot(data = ESs) +
         geom_histogram(mapping = aes(esNull, fill = ..count..), binwidth = bin_size) +
         scale_fill_gradient(name = legend_name_fill, low = "#99FF33", high = "#009900") +
         # ggplot2 Quick Reference: linetype: http://sape.inf.usi.ch/quick-reference/ggplot2/linetype
         # remove extra legends: http://stackoverflow.com/questions/11714951/remove-extra-legends-in-ggplot2
         geom_vline(mapping = aes(xintercept = ES, linetype = 'k'), colour = "#003300", size = 1) +
         geom_text(mapping = aes(x = ES, y = 0, label = paste('ES(S). Rank', ES.rank, sep =': ')), alpha = 0.4, size = 8, colour = "#003300", hjust = -1, vjust = 1.5, angle = 90) +
         scale_linetype_manual(legend_name_linetype, values = c('k' = 2), labels = c(round(ES, 2))) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_boxplot_varianceRank_boolenhancer <- function(data) {
    Title <- 'Boxplot: -log10(p-value) Rank'
    Xlab  <- "boolean enhancer"
    Ylab  <- "-log10(p-value) Rank"
    legend_name_fill <- "boolean enhancer"
    data[, rIndex := data[,.I]]
    g <- ggplot(data = data) +
         geom_boxplot(mapping = aes(boolenhancer, rIndex, fill = boolenhancer)) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

plot_pHitMiss <- function(data) {
    Title <- 'Scatterplot: pHits/pMiss v.s. -log10(p-value) Rank'
    Xlab  <- "-log10(p-value) Rank"
    Ylab  <- "pHitMiss"
    data[, rIndex := data[,.I]]
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = rIndex, y = pHitMiss, colour = boolenhancer), size = 1) +
         labs(title = Title, x = Xlab, y = Ylab) +
         get_Theme()
    return(g)
}

draw_plots12 <- function(g1, g2, outfile, Width = 2400, Height = 2112) {
    Text <- paste(gender, pheno, chrN, db, sep = ', ')
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(1, 5, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 40), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(g1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(g2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    dev.off()
}

draw_plots12345678 <- function(g1, g2, g3, g4, g5, g6, g7, g8, outfile, Width = 4500, Height = 4032) {
    Text <- paste(gender, pheno, chrN, db, sep = ', ')
    # plot distbution ratio here is: Width : Height = (2400/2*3) : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(5, 5, heights = unit(c(1, 5, 5, 5, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 100), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4))
    print(g1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:4))
    print(g2, vp = viewport(layout.pos.row = 2, layout.pos.col = 5))
    print(g3, vp = viewport(layout.pos.row = 3, layout.pos.col = 5))
    print(g4, vp = viewport(layout.pos.row = 4, layout.pos.col = 5))
    print(g5, vp = viewport(layout.pos.row = 5, layout.pos.col = 5))
    print(g6, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:4))
    print(g7, vp = viewport(layout.pos.row = 4, layout.pos.col = 1:4))
    print(g8, vp = viewport(layout.pos.row = 5, layout.pos.col = 1:4))
    dev.off()
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
}

## run ##
run <- function() {
    data1 <- read_file(infile1)  # GSEA cal result
    data2 <- read_file(infile2)  # esNulls 
    data3 <- read_file(infile4)  # EsRank
    loopList <- getLoopOverList(data1)
#    esNulls <- read_file(infile2, f22(f21(infile2, '[.]', 1), '/', 1))
#    ES.rank <- read_file(infile3, f22(f21(infile3, '[.]', 1), '/', 1))$rankES
#    ES <- get_ES(data)
#    g1 <- plot_running_sum(data)
#    g2 <- plot_histogram_variance_boolenhancer(data)
#    g3 <- plot_histogram_esNull_ES(esNulls, ES, ES.rank)
    #draw_plots12(g1, g3, outfile1)   
#    g4 <- plot_boxplot_varianceRank_boolenhancer(data)
#    g5 <- plot_effect_variance(data)
#    g6 <- plot_pHitMiss(data)
#    draw_plots12345678(g1, g2, g3, g4, g5, g6, g7, g8, outfile2)
}

run()



