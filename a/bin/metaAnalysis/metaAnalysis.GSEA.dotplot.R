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
infile1  <- args[1]  # GSEA results 
outfile1 <- args[2]


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
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

plot_runningSum <- function(data, pheno, db, lowFill, highFill, lineColor) {
    Title <- paste(pheno, db, sep = ", ")
    p <- ggplot(data, aes(x = .I, ymin = pHitMiss, ymax = pHitMiss))  + 
         theme_dose() + 
         xlab("Position in the Ranked List of Genes")
    if (by == "runningScore" || by == "all") {
        p.res <- p + geom_linerange(colour = "#DAB546")
        p.res <- p.res + geom_line(aes(y = accSum))
        es.dt <- data[boolenhancer == TRUE]
#        p.res <- p.res + geom_vline(data = es.dt, aes(xintercept = accSum), colour = "#FA5860", linetype = "dashed")
        p.res <- p.res + ylab("Runing Enrichment Score")
        p.res <- p.res + geom_hline(aes(yintercept = 0))
    }
    if (by == "position" || by == "all") {
        df2 <- data.frame(pos = which(p$data$position == 1))
        p.pos <- p + geom_vline(data = df2, aes(xintercept = pos),
             colour = "#DAB546", alpha = 0.3)
        p.pos <- p.pos + geom_line(aes(y = geneList), colour = "red")
        p.pos <- p.pos + ylab("Phenotype")
        p.pos <- p.pos + geom_hline(aes(yintercept = 0))
    }
    if (by == "runningScore")
        return(p.res)
    if (by == "position")
        return(p.pos)
    p.pos <- p.pos + xlab("") + theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
    p.res <- p.res + theme(axis.title.x = element_text(vjust = -0.3))
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(p.pos, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p.res, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    invisible(list(runningScore = p.res, position = p.pos))
}


Arguments
gseaResult gseaResult object
geneSetID geneSet ID
by one of "runningScore" or "position"
title plot title

    bin_size <- 1
    Xlab  <- paste("NbSigSnp.  bin size:", bin_size)
    Ylab  <- "Frequency"
    legend_name_fill     <- "NbSigSnps(Permutaion)"
    legend_name_linetype <- "NbSigSnp(Original)"
    xInterceptLabel      <- paste('NbSigSnp(Original)', NbSigSnpOriginal, sep =': ')
    g <- ggplot(data = NbSigSnpsPermutaion) +
         geom_histogram(mapping = aes(NbSigSnp, fill = ..count..), binwidth = bin_size) +
         scale_fill_gradient(name = legend_name_fill, low = lowFill, high = highFill) +
         geom_vline(mapping = aes(xintercept = NbSigSnpOriginal, linetype = 'a'), colour = lineColor, size = 10) +    
#         geom_text(mapping = aes(x = NbSigSnpOriginal, y = 0, label = xInterceptLabel), alpha = 0.4, size = 8, colour = lineColor, hjust = -1, vjust = 1.5, angle = 90) +
         scale_linetype_manual(legend_name_linetype, values = c('a' = 2), labels = NbSigSnpOriginal) +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 75),
               axis.title.x = element_text(size = 45),
               axis.title.y = element_text(size = 45),
               strip.text.x = element_text(size = 35),
               legend.text  = element_text(size = 35),
               legend.title = element_text(size = 35),
               axis.text.x  = element_text(size = 35),
               axis.text.y  = element_text(size = 35),
               legend.position = "none")
    return(g)
}

run <- function() {
    data <- read_file(infile1) 
    Data <- as(as(data, 'data.frame'), 'DataFrame')
}

run()



