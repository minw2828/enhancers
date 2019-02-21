# /usr/bin/rscript

## Description:
## This script finds overlap between SNPs
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 25 November 2016
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name>  <infile1> <outfile1>


# get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
outfile1 <- args[2] 


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables #
# <- f22(outfile1, '[.]', 4)[1]


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

venn_diagram <- function(data) {
    Sources <- unique(data$Source)
    n <- length(Sources)
    catColours <- colorRampPalette(brewer.pal(8, "Set1"))(n)
    Title    <- 'All SNPs: Villar (H3K27ac & H3K4me3) \n& Zhao H3K4me3 (tender & tough)'
    Subtitle <- ''
#paste(Sources, collapse = ', ')
    venn <- venn.diagram(x = lapply(Sources, function(x) data[Source == x]$snpName), 
                         filename = NULL, height = 3000, width = 3000, resolution = 500,
                         imagetype = "png", units = "px", compression = "lzw", na = "stop",
                         main = Title, sub = Subtitle,
                         main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif",
                         main.col = "black", main.cex = 5, main.just = c(0.5, 1),
                         sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif",
#                         sub.col = "black", sub.cex = 2, sub.just = c(0.5, 1),
                         sub.cex = 0,
                         category.names = Sources, cat.col = catColours,
                         cat.pos = c(-50, -40, 40, 50), cat.just = list(c(0, 0), c(-0.2, 0), c(1, 0), c(0.7, 0)), cat.cex = rep(2.5, n),
                         force.unique = TRUE, print.mode = "raw", sigdigs = 3,
                         direct.area = FALSE, area.vector = 0, hyper.test = FALSE,
                         total.population = NULL, lower.tail = TRUE,
                         col = catColours[c(1,3,4,2)], fill = catColours, alpha = rep(0.4, n),
                         cex = 4, scaled = TRUE)
    return(venn)
}

draw_plots <- function(g, outfile) {
    nr <- 1; nc <- 1
    Width  <- 1200 * nc
    Height <- 960 * nr
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(rep(5, nr)), "null"))))
    grid.draw(g)
    popViewport()
    dev.off()
}



run <- function() {
    data1 <- read_file(infile1)
    g <- venn_diagram(data1)
    draw_plots(g, outfile1)
}

run()

