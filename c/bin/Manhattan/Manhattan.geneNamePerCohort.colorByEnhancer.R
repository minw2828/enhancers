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
infile1  <- args[1]  # meta-analysis result within gene region 
infile2  <- args[2]  # enhancer sites within gene region 
outfile1 <- args[3]
# infile1 <- '/tmp//chr_pos_pheno_pmeta.tab'
# infile2 <- '/tmp//chr_pos_source.2016-01-19.a.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/Manhattan.geneName.colorByEnhancer/output.DGAT1.10.3.b.2016-01-19.a.png'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
geneName     <- f22(outfile1, '[.]', 7)[1]
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

getLoopOverList <- function(data1, data2) {
    phenos  <- unique(data1$pheno)
    chrNs   <- unique(data1$chr)
    dbs     <- unique(data2$source)
    dt <- setDT(expand.grid(phenos, chrNs, dbs))
    setnames(dt, c('pheno', 'chrN', 'db'))
    return(dt[])
}

mergeDts <- function(data1, data2) {
    setkey(data1, chr,pos)
    setkey(data2, chr,pos)
    mdt <- merge(data1, data2, all.x = TRUE)
    mdt[!is.na(source), boolEnhancer := TRUE]
    mdt[is.na(source), boolEnhancer := FALSE]
    mdt[, source := NULL]
    return(mdt)
}

getManhattanPos <- function(data) {
    n <- nchar(as.character(max(data$pos)))
    data[, manhattanPos := paste(sprintf("%02d", chr), sprintf(paste("%", n, "d", sep = ''), pos), sep = ':')]
    return(data)
}

log_transform <- function(pmeta, logBase) {
    return(-log(pmeta, logBase))
}

plot_Manhattan_geneName_db <- function(data, Pheno, chrN, Db, base, verysignifline, genomewideline, suggestiveline) {
    Title <- paste(Pheno, Db, sep =', ')
    Xlab  <- paste('chromosome', chrN, sep = ': ')
    Ylab  <- '-log10(p-value)'
    legend_name_colour <- "p-value thresholds & boolean enhancer: "
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = manhattanPos, y = log_transform(pmeta, base), colour = boolEnhancer), size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(verysignifline, base), colour = as.factor(verysignifline)), linetype = 4, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(genomewideline, base), colour = as.factor(genomewideline)), linetype = 2, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(suggestiveline, base), colour = as.factor(suggestiveline)), linetype = 1, size = 0.5) +
         scale_x_discrete(breaks = as.numeric(levels(as.factor(data$chr))), labels = as.numeric(levels(as.factor(data$chr)))) +
         scale_alpha(guide = "none") +
         guides(colour = guide_legend(title = legend_name_colour, nrow = 1)) +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 40),
               axis.title.x = element_text(size = 35),
               axis.title.y = element_text(size = 35),
               legend.text  = element_text(size = 30),
               legend.title = element_text(size = 30),
               legend.position = "bottom",
               axis.text.x  = element_text(size = 30),
               axis.text.y  = element_text(size = 30))
    return(g)
}

get1plot <- function(data1, data2, loopList, index, logBase, verysignifline, genomewideline, suggestiveline) {
    oneRow <- loopList[index]
    dat1 <- data1[pheno == oneRow$pheno]
    dat2 <- data2[source == oneRow$db]
    data <- mergeDts(dat1, dat2)
    getManhattanPos(data)
    g <- plot_Manhattan_geneName_db(data, oneRow$pheno, oneRow$chrN, oneRow$db, logBase, verysignifline, genomewideline, suggestiveline)
    return(g)
}

draw_plots <- function(data1, data2, loopList, geneName, subPlotIndex, nc, logBase, verysignifline, genomewideline, suggestiveline, outfile) {
    Title    <- paste('Manhattan plot: Meta-analysis. Zoom in to', geneName, 'gene region')
    nr <- nrow(loopList)/(nc-1) + 1
    Width  <- 3000 * ((nc-1) + 1/40*2)
    Height <- 750 * ((nr-1) + 1/10*2)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(2, rep(10, nr-1)), rep("inches", nr)), widths = unit(c(2, rep(40, nc-1)), rep("inches", nc)))))
    grid.text(Title, gp = gpar(fontsize = 55, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = nc))
    grid.text(subPlotIndex, gp = gpar(fontsize = 40, fontface = "bold"), vp = viewport(layout.pos.row = 1:nr, layout.pos.col = 1))
    for (i in seq(2, nr)) {
        index <- i - 1
        g <- get1plot(data1, data2, loopList, index, logBase, verysignifline, genomewideline, suggestiveline)
        print(g, vp = viewport(layout.pos.row = i, layout.pos.col = nc))
    }
    dev.off()
}

run <- function() {
    data1 <- read_file(infile1) # GWAS result within gene region
    data2 <- read_file(infile2) # enhancer sites within gene region 
    loopList <- getLoopOverList(data1, data2)
    draw_plots(data1, data2, loopList, geneName, subPlotIndex, nc, logBase, verysignifline, genomewideline, suggestiveline, outfile1)
}

run()


