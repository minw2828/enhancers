# /usr/bin/rscript

## Description:
## This script plots a venn diagram of significant variants from meta-analysis results across 
## all database cohorts.
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
## Rscript <module_name> <meta-analysis> <enhancer snp> <venn diagram>


# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3] # venn diagram of significant variants per enhancer set
# infile1  <- '/tmp//snpName_pheno_pvalue.10e-08.2016-01-19.c.csv'
# infile2  <- '/tmp//snpName_db.csv'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/vennDiagram.metaAnalysis.sigSnpWithinEnhancer/vennDiagram.metaAnalysis.sigSnpWithinEnhancer.2016-01-19.c.png'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(infile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(infile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
outpath  <- file.path(path_pre, 'analyses', date, psf, 'out')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pmetaSig <- f22(outfile1, '[.]', 8)[1]
mainCex    <- as.numeric(f22(outfile1, '[.]', 7)[1])
subMainCex <- as.numeric(f22(outfile1, '[.]', 6)[1])
Weight     <- as.numeric(f22(outfile1, '[.]', 5)[1])
Height     <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analyses #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getBoolenhancerPerDb <- function(data1, data2) {
    data1[snpName %in% data2$snpName, boolEnhancer := TRUE]
    data1[!(snpName %in% data2$snpName), boolEnhancer := FALSE]
    return(data1)
}

getSigSnpPerDb <- function(data1, data2, pmetaSig, Db) {
    pmetaSig <- as.numeric(pmetaSig)
    dat1 <- data1[pmeta <= pmetaSig]
    dat2 <- data2[db == Db]
    getBoolenhancerPerDb(dat1, dat2)
    res <- dat1[boolEnhancer == TRUE]
    return(res)
}

venn_diagram <- function(data1, data2, pmetaSig, mainCex, subMainCex, Weight, Height) {
    Title    <- paste('Meta-analysis significant putative bovine enhancer variants')
    Subtitle <- paste('significance level', pmetaSig, sep = ': ' )
    dbs <- unique(data2$db)
    n <- length(dbs)
    colours <- colorRampPalette(brewer.pal(n, "Set1"))(n)
    X <- lapply(dbs, function(x) getSigSnpPerDb(data1, data2, pmetaSig, x)$snpName)
    names(X) <- dbs
    venn <- venn.diagram(x = X, 
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
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    g <- venn_diagram(data1, data2, pmetaSig, mainCex, subMainCex, Weight, Height)
    draw_plots(g, Weight, Height, outfile1)
}

run()


