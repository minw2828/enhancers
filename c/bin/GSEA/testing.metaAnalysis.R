# /usr/bin/rscript

## Description:
## This script is used to test if another script is performed correctly.
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 01 July 2016
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name>


scriptName <- 'metaAnalysis.R'

# Before perceeding 
# comment run() in the script being tested 

# infile1 <- '/tmp//gender_pheno_snpName_effect_pvalue.14.2016-01-19.c.tab'
# infile2 <- '/tmp//gender_pheno_snpName_freqAllele0_freqAllele1_freqAlleleMinor_freqGeno0_freqGeno1_freqGeno2.14.2016-01-19.c.tab'
# infile3 <- '/group/dairy/Min/geno2pheno/data/phenotype/Iona_20160115/pheno_idInt_id_sex_FY_MY_PY_EDC_YoB_Breed.txt'
# nCows    <- as.numeric('12395')
# nBulls   <- as.numeric('4186')
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis/scatterplots.distributionOfZvalues.14.0.1e-01.1e-02.2016-01-19.c.png'

# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, scriptName))

# analysis #
getMeanSdPhenosByGender <- function(data) {
    # The third ('sex') is 1 for cows and 2 for bulls.
    # The next three are for Fat, Milk and Prot Yield in that order.
    data[sex == 1, gender := 'cow']
    data[sex == 2, gender := 'bull']
    dat <- data[, .(gender, FY, MY, PY)]
    melted <- melt(dat, id.vars = "gender")
    DT <- setDT(ddply(melted, .(gender, variable), summarize, meanPheno = mean(value), sdPheno = sd(value)))
    setnames(DT, 'variable', 'pheno')
    return(DT)
}

combineGenoPheno <- function(data1, data2) {
    setkey(data1, gender,pheno)
    setkey(data2, gender,pheno)
    mdt <- merge(data1, data2, all.x = TRUE)
    return(mdt)
}

getLabsTheme <- function(Title, Xlab, Ylab) {
    return(labs(title = Title, x = Xlab, y = Ylab) + 
           theme(plot.title = element_text(face = "bold", size = 25),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20),
                 strip.text.x = element_text(size = 20),
                 axis.text.x= element_text(size = 20),
                 axis.text.y = element_text(size = 20),
                 legend.position = "none"))
}

histogram_col_acrossPheno <- function(data, col) {
    bin_size <- 10000
    Title <- paste('Histogram:', col, 'across cohorts on chromsome', chrN)
    Xlab  <- col
    Ylab  <- 'frequency'
    g <- ggplot(data = data) +
         geom_histogram(mapping = aes(get(col), fill = pheno), binwidth = bin_size) +
         facet_wrap(~ pheno, ncol = 3) + getLabsTheme(Title, Xlab, Ylab)
    return(g)
}

scatterplot_xValues_yValues <- function(data, xValues, yValues) {
    Title <- paste('Scatterplot:', yValues, 'over', xValues)
    Xlab  <- xValues
    Ylab  <- yValues
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = get(xValues), y = get(yValues), colour = pheno), size = 1) +
         facet_wrap(~ pheno, ncol = 3) + getLabsTheme(Title, Xlab, Ylab)
    return(g)
}

scatterplot_variance_across <- function(DATA, SnpNames) {
    data <- DATA[snpName %in% SnpNames]
    data[, cohort := paste(gender, pheno, sep = '-')]
    Title <- 'Scatterplot: variance per SNP on across 6 cohorts'
    Xlab  <- 'cohorts'
    Ylab  <- 'variance'
    legend_name_fill <- "enhancer boolean"
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = cohort, y = variance, colour = snpName), size = 10) +
         facet_wrap(~ snpName, ncol = 3) + getLabsTheme(Title, Xlab, Ylab)
    return(g)
}

draw_plots12 <- function(g1, g2, outfile, Width = 2400, Height = 1920) {
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(5, 5), "null"))))
    print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(g2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    dev.off()
}

scatterplot_twoPQ_across <- function(DATA, SnpNames) {
    data <- DATA[snpName %in% SnpNames]
    data[, cohort := paste(gender, pheno, sep = '-')]
    Title <- 'Scatterplot: variance per SNP on across 6 cohorts'
    Xlab  <- 'cohorts'
    Ylab  <- 'variance'
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = cohort, y = twoPQ, colour = snpName), size = 10) +
         facet_wrap(~ snpName, ncol = 3) + getLabsTheme(Title, Xlab, Ylab)
    return(g)
}

scatterplot_SE_across <- function(DATA, snpNames) {
    data <- DATA[snpName %in% snpNames]
    Title <- 'Scatterplot: variance per SNP on across 6 cohorts'
    Xlab  <- 'cohorts'
    Ylab  <- 'variance'
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = cohort, y = SE, colour = snpName), size = 10) +
         facet_wrap(~ snpName, ncol = 3) + getLabsTheme(Title, Xlab, Ylab)
    return(g)
}

histogram_se_acrossCohort <- function(data) {
    data[, cohort := paste(gender, pheno, sep = '-')]
    bin_size <- 0.0001
    Title <- paste('Histogram: SE across cohorts on chromsome', chrN)
    Xlab  <- 'se'
    Ylab  <- 'frequency'
    g <- ggplot(data = data) +
         geom_histogram(mapping = aes(se, fill = cohort), binwidth = bin_size) +
         facet_wrap(~ cohort, ncol = 3) + getLabsTheme(Title, Xlab, Ylab)
    return(g)
}

getTopNMaxVariance <- function(data, topN) {
    L <- data[data[, .I[var %in% head(sort(var, decreasing = TRUE), topN)]]][order(var, decreasing = TRUE)]
    return(L)
}

getTopNMaxVarianceByCohort <- function(data, topN) {
    L <- data[, .SD[variance %in% head(sort(variance, decreasing = TRUE), topN)], by = cohort]
    return(L)
}

overlapsSnpNameAcrossCohorts <- function(L, LC) {
    DT2 <- setDT(ddply(LC, .(snpName), nrow))
    snpNames <- DT2[DT2$V1 == 6]$snpName 
    L_snpNames <- L[snpName %in% snpNames]$snpName
    LC[snpName %in% L_snpNames][order(variance, decreasing = TRUE)]$snpName
}

venndiagram_meta_ori <- function(data1, data2) {
#    data1 <- data_ori[pvalue <= 10e-8]$snpName
#    data2 <- data_meta[pmeta <= 10e-8]$snpName
    venn.plot <- draw.pairwise.venn(
        area1 <- length(unique(data1)),
        area2 <- length(unique(data2)),
        cross.area <- length(intersect(data1, data2)),
        category = c("GWAS", "meta-analysis"),
        fill = c("blue", "red"),
        lty = "blank",
        cex = 2,
        cat.cex = 2
     #   cat.pos = c(285, 105),
     #   cat.dist = 0.09,
     #   cat.just = list(c(-1, -1), c(1, 1)),
     #   ext.pos = 30,
     #   ext.dist = -0.05,
     #   ext.length = 0.85,
     #   ext.line.lwd = 2,
     #   ext.line.lty = "dashed"
    )
    return(venn.plot)
}

draw_plots12 <- function(venn.plot, outfile, n, Width = 1200 * n, Height = 960 * n) {
    Text <- paste('Venn Diagram: Overlaps of significant variants in GWAS and meta-analysis.')
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 30), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(grid.draw(venn.plot), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    dev.off()
}

run <- function() {
    # get all original inputs 
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    data <- calculate_varianceOfSnpEffect(data1, data2)
    loopList <- getLoopOverList(vts)
    # get plot color 
    Fill <- loopList[1]$Fill
    Colour <- loopList[1]$Colour
    # Add phenotype values 
    cns <- unlist(strsplit(f21(f22(infile3, '/', 1), '[.]', 1), '_'))
    data3 <- setnames(fread(infile3), c('phenoId', 'intId', cns[4:length(cns)]))
    data3.summary <- getMeanSdPhenosByGender(data3)
    Data <- combineGenoPheno(data, data3.summary)
    Data[, cohort := paste(gender, pheno, sep = '-')]
    # filter by NA and vt 
    vt = 0.01
    dat <- filterCrossAllCohorts(Data, vt)
    getStandardErrors(dat, nCows, nBulls)
    dt <- getMetalStatisticZ(dat)
    getVarianceFromMetalStatistical(dt)
    # get top N maximum variances 
    topN <- 20 
    L <- getTopNMaxVariance(dt, topN)
    # get top N maximum variance by cohort
    LC <- getTopNMaxVarianceByCohort(dat, topN)

    g <- scatterplot_variance_across(Data, L$snpName)
    outfile <- paste('top', topN, 'VariancesAcrossCohorts.png', sep = '')
    f131(outfile, g, 1200*3, 960*3)

    g <- scatterplot_zvalues(dt, vt, Fill, Colour)
    f13('zvaluesOverVariance.png', g)

    data[, twoPQ := 2 * freqAllele0 * freqAllele1]
    g <- scatterplot_twoPQ_across(Data, L$snpName)
    outfile <- paste('top', topN, 'twoPQacrossCohorts.png', sep = '')
    f131(outfile, g, 1200*3, 960*3)

    g <- scatterplot_SE_across(dt, L$snpName)
    outfile <- paste('top', topN, 'seAcrossCohorts.png', sep = '')
    f131(outfile, g, 1200*3, 960*3)

    g <- histogram_col_acrossCohort(dat, 'se')
    outfile <- 'seAcrossCohorts.png'
    f131(outfile, g, 1200, 960)

    target <- dat[variance %in% head(sort(variance, decreasing = T), 1000)]$snpName
    DAT <- dat[snpName %in% target]
    g <- histogram_se_acrossCohort(DAT)
    outfile <- 'seAcrossCohorts.png'
    f131(outfile, g, 1200*3, 960*2)

    # plot histogram of se  
    g <- histogram_col_acrossPheno(dat1, 'se')
    f131('sePerPheno.png', g, 1200*3, 960)

    # plot zvalues over se 
    g <- scatterplot_xValues_yValues(dat1, 'varmeta', 'zmeta')
    f131('seZPerPheno.png', g, 1200*3, 960)

    # plot se over effect 
    g <- scatterplot_xValues_yValues(dat1, 'effect', 'se')
    f131('effectSePerPheno.png', g, 1200*3, 960)

    data1 <- data_ori[pvalue <= 10e-8]$snpName
    data2 <- data_meta[pmeta <= 10e-8]$snpName
    draw_plots12(venndiagram_meta_ori(data1, data2), 'venndiagram.ori_meta.png', 1)
}

run()

