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
mainCex    <- as.numeric(f22(outfile1, '[.]', 7)[1])
subMainCex <- as.numeric(f22(outfile1, '[.]', 6)[1])
Weight     <- as.numeric(f22(outfile1, '[.]', 5)[1])
Height     <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

transformData <- function(data, dbs) {
    res <- lapply(dbs, function(x) data[db == x]$snpName)
    names(res) <- dbs
    return(res)
}

venn_diagram <- function(data, dbs, mainCex, subMainCex, Weight, Height) {
    Title <- 'Overlapping genomic intervals between every two enhancer sets'
    Subtitle <- '(measure by total base pairs)'
    n <- length(dbs)
    colours <- colorRampPalette(brewer.pal(n, "Set1"))(n)
    venn <- venn.diagram(x = data, 
                         filename = NULL, height = Height, width = Height, resolution = 500,
                         imagetype = "png", units = "px", compression = "lzw", na = "stop",
                         main = Title, main.pos = c(0.5, 1.05), main.fontface = "plain", 
                         main.fontfamily = "serif", main.col = "black", main.cex = mainCex, 
                         main.just = c(0.5, 1),
                         sub = Subtitle, sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                         sub.fontfamily = "serif", sub.col = "grey10", sub.cex = subMainCex,
                         sub.just = c(0.5, 1),
                         category.names = dbs, cat.col = colours, cat.dist = rep(0.1, n), 
                         cat.cex = rep(7, n),
#                         cat.pos = c(-45, 135), 
                         force.unique = TRUE, print.mode = "raw", sigdigs = 3,
                         direct.area = FALSE, area.vector = 0, hyper.test = FALSE,
                         total.population = NULL, lower.tail = TRUE,
                         col = colours, fill = colours, alpha = rep(0.4, n),
                         cex = 4.5, scaled = TRUE)
    return(venn)
}

draw_plots <- function(g, Weight, Height, outfile) {
    nr <- 1; nc <- 1
    Width  <- Weight * nc
    Height <- Height * nr
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(rep(5, nr)), "null"))))
    grid.draw(g)
    popViewport() 
    dev.off()
}

run <- function() {
    data <- read_file(infile1)
    dbs <- unique(data$db)
    Data <- transformData(data, dbs)
    g <- venn_diagram(Data, dbs, mainCex, subMainCex, Weight, Height)
    draw_plots(g, Weight, Height, outfile1)
}

run ()

