# /usr/bin/rscript

## Description:
## This script generates a boxplot to visualize the variance of three priortised
## DGAT1 enhancer variants.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 17 June 2015
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

## Get Support libraries and cleaning data ##
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
outpath  <- file.path(path_pre, 'analyses', date, psf, 'out', 'GSEA')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

## read file and give it colname ##
read_file <- function(infile) {
    data <- fread(infile, header = FALSE)
    setnames(data, unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_')))
    return(data)
}

rawProcess <- function(data1, data2) {
    setkey(data1, RefSeq)
    setkey(data2, RefSeq)
    mdt <- merge(data1, data2, all.x = TRUE)
    res <- mdt[, .(Name, start, end, geneName)]
    setnames(res, 'Name', 'chr')
    return(res)
}

nonRedundant <- function(data) {
    dat1 <- data[, min(start), by = .(chr, geneName)]
    dat2 <- data[, max(end), by = .(chr, geneName)]
    setnames(dat1, 'V1', 'start')
    setnames(dat2, 'V1', 'end')
    setkey(dat1, chr, geneName)
    setkey(dat2, chr, geneName)
    mdt <- merge(dat1, dat2)
    return(mdt)
}

read_gsea_file <- function(gender, pheno, db, chr, start, end, GeneName) {
    inName <- paste('lesESnpChr_lesESnpPos_lesESnpGwasEffect_lesESnpGwasPvalue_lesESnpGseaVariance_boolenhancer_pHitMiss_accSum_annoMethod_annoType_regionStart_regionEnd_geneName_descriptionProduct', gender, pheno, sprintf("%02d", chr), db, 1000000, date, psf, 'tab', sep = '.')
    infile <- file.path(outpath, 'allNcbiGeneAroundLesESnp', gender, pheno, sprintf("%02d", chr), db, inName)
    cns  <- unlist(strsplit(f21(inName, '[.]', 1), '_'))
    data <- setnames(fread(infile), cns)
    data[, gender := gender]
    data[, pheno := pheno]
    data[, db := db]
    res <- data[geneName == GeneName][, .(gender, pheno, db, lesESnpChr, lesESnpPos, lesESnpGwasPvalue, geneName)]
    setnames(res, c('gender', 'pheno', 'db', 'chr', 'pos', 'pvalue', 'geneName'))
    return(res)
}

getPvaluesPerRecord <- function(record) {
    # record <- data[1]
    L <- list(gender = c('bull', 'cow'), 
              pheno  = c('FY', 'MY', 'PY'), 
              db     = c('VISTA', 'FANTOM5', 'dbSUPER', 'Villar_2015'))
    dt <- as.data.table(do.call(expand.grid, L))
    n <- nrow(dt)
    dat <- rbindlist(lapply(seq(n), function(x) read_gsea_file(dt[x]$gender, dt[x]$pheno, dt[x]$db, record$chr, record$start, record$end, record$geneName)))
    dat[, cohort := paste(gender, pheno, sep = '-')]
    return(dat)
}


plot_boxplot_pvalues <- function(data) {
    Title <- paste('Boxplot: Distribution of pvalues targetting', data$geneName, 'across cohorts') 
    Xlab  <- "cohorts"
    Ylab  <- "-log10(pvalues)"
    legend_name_fill <- "cohorts"
    g <- ggplot() +
         geom_boxplot(mapping = aes(cohort, -log10(pvalue), fill = cohort), data = data) + 
         facet_wrap(~ db, ncol = 4) + 
         labs(title = Title, x = Xlab, y = Ylab) +
         guides(fill = guide_legend(title = legend_name_fill)) +
         theme(plot.title = element_text(face = "bold", size = 25),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
               axis.text.x= element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.text.x= element_text(size = 20),
               axis.text.y = element_text(size = 20))
    return(g)
}

draw1Plot2File <- function(record) {
    dat <- getPvaluesPerRecord(record)
    g <- plot_boxplot_pvalues(dat)
    outname <- paste("boxplot", record$geneName, date, psf, "png", sep = ".")
    outfile <- file.path(outpath, 'boxplotVariance.unionset', outname)
    f131(outfile, g, 4800, 960)
}

drawPlots2File <- function(data) {
    n <- nrow(data)
    for (i in seq(n)) {
        draw1Plot2File(data[i])
    }
}

draw_plots123 <- function(g1, g2, g3, outfile, Width = 2400, Height = 3840) {
    # plot distbution ratio here is: Width : Height = 2400 : (960 * 4)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(5, 5, 5), "null"))))
    print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(g2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(g3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    dev.off()
}


run <- function() {
    data1 <- read_file(infile1)
    data2 <- fread(infile2)
    data <- nonRedundant(rawProcess(data1, data2))
    drawPlots2File(data)
}

run()

