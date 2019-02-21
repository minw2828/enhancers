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
infile2  <- args[2]
outfile1 <- args[3]
# infile1  <- '/tmp//pheno_snpName_effect_pmeta.10e-08.2016-01-19.a.tab'
# infile2  <- '/tmp//snpName_source.2016-01-19.a.tab'
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
pSigLevel <- as.numeric(f22(outfile1, '[.]', 6)[1])
ntimes    <- as.numeric(f22(outfile1, '[.]', 5)[1])
nc        <- as.numeric(f22(outfile1, '[.]', 4)[1])

# analysis #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

calNbSigSnp <- function(data, sigThreshold) {
   dat <- data[pmeta <= sigThreshold]
   res <- nrow(dat)
   return(res)
}

calNbSigSnpsDb <- function(data1, data2, sigThreshold) {
    setkey(data1, snpName)
    setkey(data2, snpName)
    mdt <- merge(data1, data2)
    res <- calNbSigSnp(mdt, sigThreshold)
    return(res)
}

## randomly draw $ndraw number of rows from data and calculate the number of significant SNPs within ##
drawSnps_calNbSigSnp <- function(data, ndraw, sigThreshold) {
    dat <- data[sample(.N, ndraw, replace = TRUE)]
    calNbSigSnp(dat, sigThreshold)
}

## calculate permutations ##
### n: The number of times of permutations ###
calculate_permutations <- function(data, ndraw, sigThreshold, n) {
    res <- as.data.table(replicate(n, drawSnps_calNbSigSnp(data, ndraw, sigThreshold), simplify = "vector"))
    setnames(res, 'NbSigSnp')
    return(res)
}

getLoopOverList <- function(data1, data2) {
    phenos  <- unique(data1$pheno)
    dbs     <- unique(data2$source)
    dt <- setDT(expand.grid(phenos, dbs))
    setnames(dt, c('pheno', 'db'))
    dt[, lineColor  := colorRampPalette(brewer.pal(8, "Set1"))(nrow(dt))]
    dt[, lowFill    := colorRampPalette(brewer.pal(8, "Set2"))(nrow(dt))]
    dt[, highFill   := colorRampPalette(brewer.pal(8, "Dark2"))(nrow(dt))]
    return(dt[])
}

## plot histogram of NbSigSnps(permutation) and NbSigSnp(original) ##
plot_histogram <- function(pheno, db, NbSigSnpsPermutaion, NbSigSnpOriginal, sigThreshold, lowFill, highFill, lineColor) {
    Title <- paste(pheno, db, sep = ", ")
    bin_size <- 1
    Xlab  <- paste("NbSigSnp.  bin size:", bin_size)
    Ylab  <- "Frequency"
    legend_name_fill     <- "NbSigSnps(Permutaion)"
    legend_name_linetype <- "NbSigSnp(Original)"
    xInterceptLabel      <- paste('NbSigSnp(Original)', NbSigSnpOriginal, sep =': ')
    g <- ggplot(data = NbSigSnpsPermutaion) +
         geom_histogram(mapping = aes(NbSigSnp, fill = ..count..), binwidth = bin_size) +
         scale_fill_gradient(name = legend_name_fill, low = lowFill, high = highFill) +
         geom_vline(mapping = aes(xintercept = NbSigSnpOriginal, linetype = 'k'), colour = lineColor, size = 1) +
         geom_text(mapping = aes(x = NbSigSnpOriginal, y = 0, label = xInterceptLabel), alpha = 0.4, size = 8, colour = lineColor, hjust = -1, vjust = 1.5, angle = 90) +
         scale_linetype_manual(legend_name_linetype, values = c('k' = 2), labels = NbSigSnpOriginal) +
         labs(title = Title, x = Xlab, y = Ylab) +
         theme(plot.title   = element_text(face = "bold", size = 35),
               axis.title.x = element_text(size = 30),
               axis.title.y = element_text(size = 30),
               strip.text.x = element_text(size = 30),
               legend.text  = element_text(size = 25),
               legend.title = element_text(size = 25),
               axis.text.x  = element_text(size = 25),
               axis.text.y  = element_text(size = 25))
    return(g)
}

get1plot <- function(data1, data2, loopList, index, pSigLevel, ntimes) {
    oneRow <- loopList[index]
    dat1 <- data1[pheno == oneRow$pheno]
    dat2 <- data2[source == oneRow$db]
    ndraw <- nrow(dat2)
    NbSigSnpOriginal    <- calNbSigSnpsDb(dat1, dat2, pSigLevel)
    NbSigSnpsPermutaion <- calculate_permutations(dat1, ndraw, pSigLevel, ntimes)
    g <- plot_histogram(oneRow$pheno, oneRow$db, NbSigSnpsPermutaion, NbSigSnpOriginal, pSigLevel, oneRow$lowFill, oneRow$highFill, oneRow$lineColor)
    return(g)
}

draw_plots <- function(data1, data2, loopList, pSigLevel, ntimes, nc, outfile) {
    Title    <- paste(ntimes, 'permutation test for enrichment of meta-analysis significant variants within enhancer sets')
    Subtitle <- paste('P-value Significance:', pSigLevel)
    nr <- nrow(loopList)/nc + 2
    Width  <- 1200 * nc
    Height <- 960 * ((nr-2) + 1/5*2)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(1, 0.5, rep(5, nr-2)), "null"))))
    grid.text(Title, gp = gpar(fontsize = 60, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:nc))
    grid.text(Subtitle, gp = gpar(fontsize = 50, fontface = "bold"), vp = viewport(layout.pos.row = 2, layout.pos.col = 1:nc))
    for (i in seq(3, nr)) {
        for (j in seq(1, nc)) {
            index <- nc * i + j - nc * 3
            g <- get1plot(data1, data2, loopList, index, pSigLevel, ntimes)
            print(g, vp = viewport(layout.pos.row = i, layout.pos.col = j))
        }
    }
    dev.off()
}

## run ##
run <- function() {
    data1 <- read_file(infile1) # all META results 
    data2 <- read_file(infile2) # enhancer variants
    loopList <- getLoopOverList(data1, data2)
    draw_plots(data1, data2, loopList, pSigLevel, ntimes, nc, outfile1)
}

run()



