# /usr/bin/rscript

## Description:
## This module plots up TADB permutation results from bostau8 and hg38.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 26 August 2016
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name> <infile1> <outfile1>

# Get input arguments #
args <- commandArgs(trailingOnly = FALSE)
index    <- 5
program  <- gsub('--file=', '', args[index-1])
infile1  <- args[index+1]    # nbSNPs 
outfile1 <- args[index+2]    # plot
outfile2 <- args[index+3]    # rank + fold change


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.bioconductor.R'))
source(file.path(binpath, 'ReadWriteFile.R'))
options(scipen = 999) # force not to use scienteific notation


# get global variables # 
Width   <- as.numeric(f22(outfile1, '[.]', 8)[1])
Height  <- as.numeric(f22(outfile1, '[.]', 7)[1])
pSig    <- as.numeric(f22(outfile1, '[.]', 6)[1])
ntimes  <- as.numeric(f22(outfile1, '[.]', 5)[1])
ORefseq <- f22(outfile1, '[.]', 4)[1]


# analysis #
rawProcess1 <- function(data1) {
    data1[, id := seq_len(.N), by = list(pheno, db, chr)]
    res <- setDT(ddply(data1, .(id, pheno, db), summarize, nbSnps = sum(nbSnp)))
    return(res)
}

getLoopOverList <- function(data) {
    dt <- unique(data[, .(pheno, db)])
    dt[, `:=` (lineColor = colorRampPalette(brewer.pal(8, "Set1"))(nrow(dt)),
               lowFill   = colorRampPalette(brewer.pal(8, "Set2"))(nrow(dt)),
               highFill  = colorRampPalette(brewer.pal(8, "Dark2"))(nrow(dt)))]
    return(dt[])
}

plot_histogram <- function(data, bin_size, ns, lowFill, highFill, lineColor) {
    Title <- paste(unique(data$pheno), unique(data$db), sep = ", ")
    g <- ggplot() +
         geom_histogram(mapping = aes(nbSnps, fill = ..count..), data = data[id != 1], binwidth = bin_size) +
         scale_fill_gradient(low = lowFill, high = highFill) +
         geom_vline(mapping = aes(xintercept = nbSnps, linetype = 'k'), data = data[id == 1], colour = "black", size = 5*ns) +
         labs(title = Title) +
         theme(plot.title      = element_text(face = "bold", size = 75*ns, hjust = 0.5),
               axis.title.x    = element_blank(),
               axis.title.y    = element_blank(),
               axis.text.x     = element_text(face = "bold", size = 40*ns, angle = 45, hjust = 1, vjust = 1),
               axis.text.y     = element_text(face = "bold", size = 40*ns, angle = 45, hjust = 1, vjust = 1),
               legend.position = "none",
               plot.margin = unit(c(1, 1, 1, 1), "cm"))
    return(g)
}

get1plot <- function(data, bin_size, ns, loopList, index) {
    oneRow <- loopList[index]
    dat <- data[pheno == oneRow$pheno & db == oneRow$db]
    g <- plot_histogram(dat, bin_size, ns, oneRow$lowFill, oneRow$highFill, oneRow$lineColor)
    return(g)
}

draw_plots <- function(data, loopList, bin_size, ns, ntimes, Width, Height, outfile) {
    nr <- length(unique(loopList$pheno))
    nc <- length(unique(loopList$db))
    rPlot <- 5
    Width  <- 1200 * nc
    Height <- 960 * nr
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(rep(rPlot, nr)), "null"))))
    for (i in seq(1, nr)) {
        for (j in seq(1, nc)) {
            index <- nc * i + j - nc * 1
            g <- get1plot(data, bin_size, ns, loopList, index)
            print(g, vp = viewport(layout.pos.row = i, layout.pos.col = j))
        }
    }
    dev.off()
}

getRanking <- function(data, loopList, index, ntimes) {
    oneRow <- loopList[index]
    dat <- data[pheno == oneRow$pheno & db == oneRow$db]
    sdat <- dat[order(nbSnps)]
    sdat[, Rank := .I]
    sdat[, rank := ifelse(id == 1 & Rank == ntimes + 1, paste('<', 1/ntimes, sep = ''), as.character(Rank/ntimes))]
    res <- unique(sdat[id == 1, .(pheno, db, rank)])
    return(res)
}

getFoldChange <- function(data, loopList, index, ntimes) {
    oneRow <- loopList[index]
    dat <- data[pheno == oneRow$pheno & db == oneRow$db]
    fc <- dat[id == 1]$nbSnps / mean(dat[id != 1]$nbSnps)
    dat[, foldChange := fc]
    res <- unique(dat[, .(pheno, db, foldChange)])
    return(res)
}

getResultRankingFoldChange <- function(data, loopList, ntimes, outfile, boolColname) {
    ns <- seq(nrow(loopList))
    tmp1 <- rbindlist(lapply(ns, function(x) getRanking(data, loopList, x, ntimes)))
    tmp2 <- rbindlist(lapply(ns, function(x) getFoldChange(data, loopList, x, ntimes)))
    setkey(tmp1, pheno,db)
    setkey(tmp2, pheno,db)
    mdt <- merge(tmp1, tmp2)
    write_file(mdt, outfile, boolColname)
}

run <- function() {
    data1 <- read_file(infile1)

    data <- rawProcess1(data1)

    loopList <- getLoopOverList(data)
 
    bin_size <- 1
    ns <- 1.9
    draw_plots(data, loopList, bin_size, ns, ntimes, Width, Height, outfile1)
    getResultRankingFoldChange(data, loopList, ntimes, outfile2, TRUE)
}

run()



