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
infile1  <- args[1] # lesESnp results (H3K4me3 & H3K27ac)
outfile1 <- args[2] # boxplots


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation


# Get global variables #
pSigLvl1 <- as.numeric(f22(outfile1, '[.]', 8)[1])
pSigLvl2 <- as.numeric(f22(outfile1, '[.]', 7)[1])
Width    <- as.numeric(f22(outfile1, '[.]', 6)[1])
Height   <- as.numeric(f22(outfile1, '[.]', 5)[1])
nc       <- as.numeric(f22(outfile1, '[.]', 4)[1])


# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

getLoopOverList <- function(data, nc) {
    if (nc == 2) { 
        odata <- data[order(pheno, db)]
    } else if (nc == 3) {
        odata <- data[order(db, pheno)]
    } else {
        odata <- data
    }
    dt <- unique(odata[, .(pheno, db)])
    dt[, dotColor  := colorRampPalette(brewer.pal(8, "Set1"))(nrow(dt))]
    return(dt[])
}

plot_scatterplot <- function(data, Pheno, Db, pSigLvl1, pSigLvl2, dotColor, nc, outfile) {
    Title <- paste(Pheno, Db, sep = ', ')
    Xlab  <- "\nposition of sorted SNPs"
    Ylab  <- "-log10(P)\n"
    legend_name_colour   <- "-log10(P): "
    legend_name_linetype <- "significance threshold: "
    g <- ggplot() +
         geom_point(mapping = aes(rank, minusLog10Pmeta), data = data, colour = dotColor, size = 2) +
         geom_hline(mapping = aes(yintercept = pSigLvl1, linetype = 'a', colour = as.factor(pSigLvl1)), size = 2, alpha = 0.5) +
         geom_hline(mapping = aes(yintercept = pSigLvl2, linetype = 'a', colour = as.factor(pSigLvl2)), size = 2, alpha = 0.5) +
         guides(colour   =  guide_legend(title = legend_name_colour, keywidth = 0.6, keyheight = 0.6, default.unit = "inch", nrow = 1),
                linetype = "none") + 
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 65),
               strip.text.x = element_text(size = 50),
               axis.title.x = element_text(size = 50), axis.title.y = element_text(size = 50),
               axis.text.x  = element_text(size = 50), axis.text.y  = element_text(size = 50),
               legend.title = element_text(size = 45, colour = 'grey15', family = "Courier"), legend.text  = element_text(size = 45, colour = 'grey15', family = "Courier"),
               legend.position = "bottom",
#               legend.position = "none", 
               plot.margin = unit(c(1, 1, 1, 1), "cm"))
    return(g)
}

draw_plots <- function(data, pSigLevel1, pSigLevel2, width, height, nc, outfile) {
    loopList <- getLoopOverList(data, nc)
    Title    <- 'GSEA core SNPs'
    nrTitle <- 1; nrPlot <- nrow(loopList)/nc
    nr <- nrTitle + nrPlot
    rTitle <- 1; rPlot <- 6
    Width  <- width * nc
    Height <- height * (nrPlot * rPlot + nrTitle * rTitle) / rPlot
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(rTitle, rep(rPlot, nrPlot)), "null"))))
    grid.text(Title, gp = gpar(fontsize = 80, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:nc))
    for (i in seq(2, nr)) {
        for (j in seq(1, nc)) {
            index <- nc * i + j - nc * 2
            dat <- data[pheno == loopList[index]$pheno & db == loopList[index]$db]
            dat[, rank := .I]
            g <- plot_scatterplot(dat, loopList[index]$pheno, loopList[index]$db, pSigLevel1, pSigLevel2, loopList[index]$dotColor, nc, outfile)
            print(g, vp = viewport(layout.pos.row = i, layout.pos.col = j))
        }
    }
    dev.off()
}

run <- function() {
    data <- read_file(infile1) 

    draw_plots(data, -log10(pSigLvl1), -log10(pSigLvl2), Width, Height, nc, outfile1)
}

run()



