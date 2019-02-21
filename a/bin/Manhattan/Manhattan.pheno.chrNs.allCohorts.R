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
infile1  <- args[1]  # meta-analysis result 
outfile1 <- args[2]
# infile1  <- '/tmp//chr_pos_pheno_pmeta.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/Manhattan.pheno.chrNs.allCohorts/output.PY.05_14_20.10.2.b.2016-01-19.c.png'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pheno        <- f22(outfile1, '[.]', 8)[1]
chrs         <- as.numeric(unlist(strsplit(f22(outfile1, '[.]', 7)[1], '_')))
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

subsetByGenderPheno <- function(data, Gender, Pheno) {
    res <- data[gender == Gender & pheno == Pheno]
    return(res)
}

getMinusLogBasePvalue <- function(data, logBase) {
    data[, minusLogBasePvalue := -log(x = pvalue, base = logBase)]
    return(data)
}

getManhattanPos <- function(data) {
    # a single phenotype data
    n <- nchar(as.character(max(data$pos)))
    data[, manhattanPos := paste(sprintf("%02d", chr), sprintf(paste("%", n, "d", sep = ''), pos), sep = ':')]
    return(data)
}

getColourList <- function(data, chrs) {
    # a single phenotype data
    dt <- data.table(chr = seq(1, 30))
    dt[, dotColor  := colorRampPalette(brewer.pal(12, "Paired"))(nrow(dt))]
    res <- dt[chr %in% chrs]
    setkey(data, chr)
    setkey(res, chr)
    mdt <- merge(data, res, all.x = TRUE)
    return(mdt)
}

getReadyForPlot <- function(data1, Gender, Pheno, logBase, chrs) {
    data <- subsetByGenderPheno(data1, Gender, Pheno)
    getMinusLogBasePvalue(data, logBase)
    getManhattanPos(data)
    res <- getColourList(data, chrs)
    return(res)
}

plot_Manhattan_allCohorts <- function(data1, Gender, Pheno, logBase, chrs, verysignifline, genomewideline, suggestiveline) {
    data <- getReadyForPlot(data1, Gender, Pheno, logBase, chrs)
    Title <- paste(Gender, Pheno, sep =', ')
    Xlab  <- paste('chromosome', paste(chrs, collapse = ', '))
    Ylab  <- '-log10(p-value)'
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = manhattanPos, y = minusLogBasePvalue), colour = data$dotColor, size = 0.5) +
         geom_hline(mapping = aes(yintercept = -log(x = verysignifline, base = logBase), colour = as.factor(verysignifline)), linetype = 4, size = 0.5) +
         geom_hline(mapping = aes(yintercept = -log(x = genomewideline, base = logBase), colour = as.factor(genomewideline)), linetype = 2, size = 0.5) +
         geom_hline(mapping = aes(yintercept = -log(x = suggestiveline, base = logBase), colour = as.factor(suggestiveline)), linetype = 1, size = 0.5) +
         scale_x_discrete(breaks = as.numeric(levels(as.factor(data$chr))), labels = as.numeric(levels(as.factor(data$chr)))) +
         scale_alpha(guide = "none") +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 40),
               axis.title.x = element_text(size = 35),
               axis.title.y = element_text(size = 35),
               legend.text  = element_text(size = 30),
               legend.title = element_text(size = 30),
               legend.position = "none",
               axis.text.x  = element_text(size = 30),
               axis.text.y  = element_text(size = 30))
    return(g)
}

getLoopOverList <- function(data1) {
    genders <- unique(data1$gender)
    phenos  <- unique(data1$pheno)
    dt <- setDT(expand.grid(genders, phenos))
    setnames(dt, c('gender', 'pheno'))
    return(dt[])
}

draw_plots <- function(data1, loopList, logBase, chrs, verysignifline, genomewideline, suggestiveline, nc, subPlotIndex, outfile) {
    Title    <- 'Manhattan plots: GWAS'
    Subtitle <- paste(paste('verysignifline', format(verysignifline, scientific = TRUE), sep = ':'),
                      paste('genomewideline', format(genomewideline, scientific = TRUE), sep = ':'),
                      paste('suggestiveline', format(suggestiveline, scientific = TRUE), sep = ':'), sep = '; ')
    nr <- nrow(loopList)/(nc-1) + 2
    Width  <- 7286 * ((nc-1) + 1/100*2)
    Height <- 728.6 * ((nr-2) + 1/10*3)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc + 1, heights = unit(c(2, 1, rep(10, nr-2)), rep("inches", nr)), widths = unit(c(2, rep(100, nc-1)), rep("inches", nc)))))
    grid.text(Title, gp = gpar(fontsize = 60, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:nc))
    grid.text(Subtitle, gp = gpar(fontsize = 40, fontface = "bold"), vp = viewport(layout.pos.row = 2, layout.pos.col = 1:nc))
    grid.text(subPlotIndex, gp = gpar(fontsize = 40, fontface = "bold"), vp = viewport(layout.pos.row = 1:nr, layout.pos.col = 1))
    for (i in seq(3, nr)) {
        index <- i - 2
        g <- plot_Manhattan_allCohorts(data1, loopList[index]$gender, loopList[index]$pheno, logBase, chrs, verysignifline, genomewideline, suggestiveline)
        print(g, vp = viewport(layout.pos.row = i, layout.pos.col = nc))
    }
    dev.off()
}

run <- function() {
    data1 <- read_file(infile1)
    loopList <- getLoopOverList(data1)
    draw_plots(data1, loopList, logBase, chrs, verysignifline, genomewideline, suggestiveline, nc, subPlotIndex, outfile1)
}

run()


