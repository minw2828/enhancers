# /usr/bin/rscript

## Description:
## This script performs meta-analysis by re-engineering BH's zScore.f90 script.
##   
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 06 July 2016 
## 
## Date modified and reason: 
## 13 July 2016    Use qnorm(value) to calcualte x value.
##
## Execution: 
## Rscript <module_name> 


# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variable #
pSigLevel <- as.numeric(f22(outfile1, '[.]', 5)[1])
chrN      <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

rawProcess <- function(data) {
    dat1 <- data[correction == 'before']
    dat2 <- data[correction == 'after']
    setkey(dat1, snpName,pheno)
    setkey(dat2, snpName,pheno)
    mdt <- merge(dat1, dat2)
    return(mdt)
}

getCorrelation <- function(data) {
    dt <- data.table(corEffect = cor(data$snpeffectmeta.x, data$snpeffectmeta.y), 
                     corPvalue = cor(data$pmeta.x, data$pmeta.y))
    return(dt)
}

plotCorEffect <- function(data) {
    Title <- 'Meta-analysis effects' 
    Xlab  <- 'before'
    Ylab  <- 'after'
    legend_name_fill     <- ''
    legend_name_linetype <- ''
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = snpeffectmeta.x, y = snpeffectmeta.y, colour = pheno)) + 
         geom_abline(intercept = 0, slope = 1) + 
         facet_wrap(~ pheno, nrow = 1) + 
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 45),
               axis.title.x = element_text(size = 40),
               axis.title.y = element_text(size = 40),
               strip.text.x = element_text(size = 40),
               legend.text  = element_text(size = 30),
               legend.title = element_text(size = 30),
               axis.text.x  = element_text(size = 30),
               axis.text.y  = element_text(size = 30),
               legend.position = "none")
    return(g)
}

plotCorPvalue <- function(data) {
    Title <- 'Meta-analysis p-values'
    Xlab  <- 'before'
    Ylab  <- 'after'
    legend_name_fill     <- ''
    legend_name_linetype <- ''
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = pmeta.x, y = pmeta.y, colour = pheno)) +
         geom_abline(intercept = 0, slope = 1) +
         facet_wrap(~ pheno, nrow = 1) +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 45),
               axis.title.x = element_text(size = 40),
               axis.title.y = element_text(size = 40),
               strip.text.x = element_text(size = 40),
               legend.text  = element_text(size = 30),
               legend.title = element_text(size = 30),
               axis.text.x  = element_text(size = 30),
               axis.text.y  = element_text(size = 30),
               legend.position = "none")
    return(g)
}

draw_plots <- function(data, outfile) {
    Title    <- 'Correlation before and after correcting for the effect of DGAT1 mutation'
    nrTitle <- 1; nrSubtitle <- 1; nrPlot <- 2
    nr <- nrTitle + nrSubtitle + nrPlot
    nc <- 1
    rTitle <- 1; rSubtitle <- 0.5; rPlot <- 15
    Width  <- 3600 * nc
    Height <- 960 * (nrPlot * rPlot + nrTitle * rTitle + nrSubtitle * rSubtitle) / rPlot
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(rTitle, rPlot, rSubtitle, rPlot), "null"))))
    grid.text(Title, gp = gpar(fontsize = 80, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    g <- plotCorEffect(data)
    print(g, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid.text('', gp = gpar(fontsize = 80, fontface = "bold"), vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    g <- plotCorPvalue(data)
    print(g, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
    dev.off()
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = ',',
                row.names = FALSE, col.names = TRUE)
}

run <- function() {
    data <- read_file(infile1)
    Data <- rawProcess(data)
    res <- getCorrelation(Data)
    write_file(res, outfile1)
    draw_plots(Data, outfile2)
}

run()


