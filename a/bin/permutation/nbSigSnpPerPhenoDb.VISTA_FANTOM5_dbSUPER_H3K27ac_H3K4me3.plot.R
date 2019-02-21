# /usr/bin/rscript

## Description:
## This module performs permutation test on meta-analysis results. 
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 14 July 2016
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <outfile1>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
outfile1 <- args[2]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/permutation.NbSigSnpPerPhenoDb/NbSigSnp_tpe_pheno_db.10e-08.10000.2016-01-19.a.csv'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/permutation.NbSigSnpPerPhenoDb/permutation.NbSigSnpPerPhenoDb.10e-08.10000.6.2016-01-19.a.png'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables # 
pSigLevel <- f22(outfile1, '[.]', 6)[1]
ntimes    <- as.numeric(f22(outfile1, '[.]', 5)[1])
nc        <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getLoopOverList <- function(data) {
    dt <- unique(data[, .(pheno, db)])
    dt[, lineColor  := colorRampPalette(brewer.pal(8, "Set1"))(nrow(dt))]
    dt[, lowFill    := colorRampPalette(brewer.pal(8, "Set2"))(nrow(dt))]
    dt[, highFill   := colorRampPalette(brewer.pal(8, "Dark2"))(nrow(dt))]
    return(dt[])
}

plot_histogram <- function(pheno, db, NbSigSnpsPermutaion, NbSigSnpOriginal, sigThreshold, lowFill, highFill, lineColor) {
    Title <- paste(pheno, db, sep = ", ")
    bin_size <- 1
    Xlab  <- paste("\nNb significant SNPs  (bin size:", bin_size, ')', sep = '')
    Ylab  <- "Frequency\n"
    legend_name_fill     <- "NbSigSnps(Permutaion)"
    legend_name_linetype <- "NbSigSnp(Original)"
    xInterceptLabel      <- paste('NbSigSnp(Original)', NbSigSnpOriginal, sep =': ')
    g <- ggplot(data = NbSigSnpsPermutaion) +
         geom_histogram(mapping = aes(NbSigSnp, fill = ..count..), binwidth = bin_size) +
         scale_fill_gradient(name = legend_name_fill, low = lowFill, high = highFill) +
         geom_vline(mapping = aes(xintercept = NbSigSnpOriginal, linetype = 'a'), colour = lineColor, size = 10) +    
         scale_linetype_manual(legend_name_linetype, values = c('a' = 2), labels = NbSigSnpOriginal) +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 90),
               strip.text.x = element_text(size = 75),
               axis.title.x = element_text(size = 70), axis.title.y = element_text(size = 70),
               axis.text.x  = element_text(size = 70), axis.text.y  = element_text(size = 70),
               legend.title = element_text(size = 70), legend.text  = element_text(size = 70), 
               legend.position = "none",
               plot.margin = unit(c(2, 2, 2, 2), "cm"))
    return(g)
}

get1plot <- function(data, loopList, index, pSigLevel, ntimes) {
    oneRow <- loopList[index]
    dat <- data[pheno == oneRow$pheno & db == oneRow$db]
    NbSigSnpOriginal    <- dat[tpe == 'ori']$NbSigSnp
    NbSigSnpsPermutaion <- dat[tpe != 'ori']
    g <- plot_histogram(oneRow$pheno, oneRow$db, NbSigSnpsPermutaion, NbSigSnpOriginal, pSigLevel, oneRow$lowFill, oneRow$highFill, oneRow$lineColor)
    return(g)
}

draw_plots <- function(data, loopList, pSigLevel, ntimes, nc, outfile) {
    Title    <- 'Enrichment of meta-analysis significant variants'
    Subtitle <- paste(paste(ntimes, 'repeats permutations'), paste('significance', pSigLevel, sep = ': '), sep = '; ')
    nrTitle <- 1; nrSubtitle <- 1; nrPlot <- nrow(loopList)/nc
    nr <- nrTitle + nrSubtitle + nrPlot 
    rTitle <- 1; rSubtitle <- 0.5; rPlot <- 5
    Width  <- 1500 * nc
    Height <- 1300 * (nrPlot * rPlot + nrTitle * rTitle + nrSubtitle * rSubtitle) / rPlot
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(rTitle, rSubtitle, rep(rPlot, nrPlot)), "null"))))
    grid.text(Title, gp = gpar(fontsize = 130, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:nc))
    grid.text(Subtitle, gp = gpar(fontsize = 90, fontface = "bold"), vp = viewport(layout.pos.row = 2, layout.pos.col = 1:nc))
    for (i in seq(3, nr)) {
        for (j in seq(1, nc)) {
            index <- nc * i + j - nc * 3
            g <- get1plot(data, loopList, index, pSigLevel, ntimes)
            print(g, vp = viewport(layout.pos.row = i, layout.pos.col = j))
        }
    }
    dev.off()
}

run <- function() {
    data <- read_file(infile1) 
    loopList <- getLoopOverList(data)
    draw_plots(data, loopList, pSigLevel, ntimes, nc, outfile1)
}

run()



