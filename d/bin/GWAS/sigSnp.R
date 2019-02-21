# /usr/bin/rscript

## Description:
## This script takes a csv file as input and output a plot to visualize 
## input statistics.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
##  2016 
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf>


# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]
outfile2 <- args[4] # venn diagram of significant variants per enhancer set
# infile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/sigSnp/chr_pos_effect_pvalue_gender_pheno.10e-08.2016-01-19.a.tab'

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
pvalueSig <- f22(infile1, '[.]', 4)[1]


# analyses #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(data)
}

rawProcess <- function(data, pvalueSig) {
    pvalueSig <- as.numeric(pvalueSig)
    res <- data[pvalue <= pvalueSig]
    return(res)
}

write_file <- function(content, outfile) {
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = '\t',
                row.names = FALSE, col.names = FALSE)
}

chrPos_intersectByGroup <- function(data) {
    data[, groupName := paste(gender, pheno, sep = '-')]
    data[, chr := as.numeric(gsub('Chr', '', f21(snpName, ':', 1)))]
    data[, pos := as.numeric(f21(snpName, ':', 2))]
    res <- as.data.table(Reduce(intersect, data[, .(list(unique(snpName))), groupName]$V1))
    res[, chr := as.numeric(gsub('Chr', '', f21(V1, ':', 1)))]
    res[, pos := as.numeric(f21(V1, ':', 2))]
    data[, groupName := NULL]
    setnames(res, 'V1', 'snpName')
    return(res[])
}

intersectByGroup <- function(data, res) {
    result <- data[snpName %in% res$snpName]
    return(result)
}

nbSigSnpInEnhancer <- function(data1, data2, db) {
    dat1 <- unique(data1[, .(snpName)])
    dat2 <- unique(data2[source == db, .(snpName)])
    res <- dat2[dat2$snpName %in% dat1$snpName]
    return(res)
}

venn_diagram <- function(data1, data2) {
    dt1 <- nbSigSnpInEnhancer(data1, data2, 'VISTA')
    dt2 <- nbSigSnpInEnhancer(data1, data2, 'FANTOM5')
    dt3 <- nbSigSnpInEnhancer(data1, data2, 'dbSUPER')
    dt4 <- nbSigSnpInEnhancer(data1, data2, 'Villar_2015')

    fo12 <- Reduce(intersect, list(dt1$snpName, dt2$snpName))
    fo13 <- Reduce(intersect, list(dt1$snpName, dt3$snpName))
    fo14 <- Reduce(intersect, list(dt1$snpName, dt4$snpName))
    fo23 <- Reduce(intersect, list(dt2$snpName, dt3$snpName))
    fo24 <- Reduce(intersect, list(dt2$snpName, dt4$snpName))
    fo34 <- Reduce(intersect, list(dt3$snpName, dt4$snpName))
    fo123 <- Reduce(intersect, list(dt1$snpName, dt2$snpName, dt3$snpName))
    fo124 <- Reduce(intersect, list(dt1$snpName, dt2$snpName, dt4$snpName))
    fo134 <- Reduce(intersect, list(dt1$snpName, dt3$snpName, dt4$snpName))
    fo234 <- Reduce(intersect, list(dt2$snpName, dt3$snpName, dt4$snpName))
    fo1234 <- Reduce(intersect, list(dt1$snpName, dt2$snpName, dt3$snpName, dt4$snpName))

    colours <- c("maroon1", "gold1", "deepskyblue1", "darkseagreen1")
    venn.plot <- draw.quad.venn(area1 = nrow(dt1),
                                area2 = nrow(dt2),
                                area3 = nrow(dt3),
                                area4 = nrow(dt4),
                                n12 = length(fo12),
                                n13 = length(fo13),
                                n14 = length(fo14),
                                n23 = length(fo23),
                                n24 = length(fo24),
                                n34 = length(fo34),
                                n123 = length(fo123),
                                n124 = length(fo124),
                                n134 = length(fo134),
                                n234 = length(fo234),
                                n1234 = length(fo1234),
                                category = c("VISTA", "FANTOM5", "dbSUPER", "Villar 2015"),
                                col = c("maroon1", "deepskyblue1", "darkseagreen1", "gold1"), # Must adjust colour mannual. The auto version of wrong, misleading, because colours do not match annotation.
                                fill = colours,
                                alpha = rep(0.5, 4),
#                                cex = rep(2, 15),
                                cex = rep(6, 15),
#                                cat.cex = rep(2, 4),
                                cat.cex = rep(5, 4),
                                cat.col = colours)
    return(venn.plot)
}

draw_plots12 <- function(venn.plot, pvalueSig, outfile, Width = 1200, Height = 960) {
    Title    <- paste('GWAS significant putative bovine enhancer variants')
    Subtitle <- paste('significance level', pvalueSig, sep = ': ')
    
    nr <- 4; nc <- 1
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr, nc, heights = unit(c(1, 1, 1, 5), "null"))))
    grid.text(Title, gp = gpar(fontsize = 50), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.text(Subtitle, gp = gpar(fontsize = 40), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid.text(" ", gp = gpar(fontsize = 25), vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(grid.draw(venn.plot), vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
    dev.off()
}

venn_diagram <- function(data1, data2, pvalueSig) {
    Title    <- paste('GWAS significant putative bovine enhancer variants')
    Subtitle <- paste('significance level', pvalueSig, sep = ': ' )
    venn <- venn.diagram(x = list(VISTA   = nbSigSnpInEnhancer(data1, data2, 'VISTA')$snpName,
                                  FANTOM5 = nbSigSnpInEnhancer(data1, data2, 'FANTOM5')$snpName,
                                  dbSUPER = nbSigSnpInEnhancer(data1, data2, 'dbSUPER')$snpName,
                                  Villar  = nbSigSnpInEnhancer(data1, data2, 'Villar_2015')$snpName),
                         filename = NULL, height = 3000, width = 3000, resolution = 500,
                         imagetype = "png", units = "px", compression = "lzw", na = "stop",
                         main = Title, sub = Subtitle,
                         main.pos = c(0.5, 1.05), main.fontface = "plain", main.fontfamily = "serif",
                         main.col = "black", main.cex = 4, main.just = c(0.5, 1),
                         sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif",
                         sub.col = "black", sub.cex = 4, sub.just = c(0.5, 1),
                         category.names = c("VISTA", "FANTOM5", "dbSUPER", "Villar"), cat.col = c("maroon", "goldenrod", "deepskyblue4", "seagreen4"),
                         cat.pos = c(-50, -40, 40, 50), #cat.dist = rep(-0.05, 4),
                         cat.just = list(c(0, 0), c(-0.5, 0), c(1, 0), c(0.7, 0)), cat.cex = rep(5, 4),
                         force.unique = TRUE, print.mode = "raw", sigdigs = 3,
                         direct.area = FALSE, area.vector = 0, hyper.test = FALSE,
                         total.population = NULL, lower.tail = TRUE,
                         col = c("maroon1", "deepskyblue1", "seagreen3", "gold1"), fill = c("maroon1", "gold1", "deepskyblue1", "seagreen3"), alpha = 0.4,
                         cex = 6, scaled = TRUE)
    return(venn)
}

draw_plots <- function(g, outfile) {
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    nr <- 1;     nc <- 1
    rhText <- 1; rhPlot <- 5
    nText <- 0;  nPlot <- 1
    Width  <- 1500 * nc
    Height <- 960 * (rhText * rhText + rhPlot * nPlot)/rhPlot
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    plot.new() # start new page
    gl <- grid.layout(nrow = nr, ncol = nc, heights = unit(c(rep(rhText, nText), rep(rhPlot, nPlot)), "null")) # setup layout
    # grid.show.layout(gl)
    vp.1 <- viewport(layout.pos.row = 1, layout.pos.col = 1) # setup viewports
    pushViewport(viewport(layout = gl)) # init layout
    pushViewport(vp.1) # access the first position
    grid.draw(g)
    popViewport() # done with the second viewport
    dev.off()
}


run <- function() {
    data <- read_file(infile1)
    data1 <- rawProcess(data, pvalueSig)

#    nrow(data1[gender == 'bull' && pheno == 'FY'])
#    nrow(data1[gender == 'bull' && pheno == 'MY'])
#    nrow(data1[gender == 'bull' && pheno == 'PY'])
#    nrow(data1[gender == 'cow' && pheno == 'FY'])
#    nrow(data1[gender == 'cow' && pheno == 'MY'])
#    nrow(data1[gender == 'cow' && pheno == 'PY'])

    res <- chrPos_intersectByGroup(data1)
    resIntersection <- intersectByGroup(data1, res)
    write_file(resIntersection, outfile1)

    data2 <- read_file(infile2)
    data2[, snpName := paste('Chr', chr, ':', pos, sep = '')]

    nbSigSnpInEnhancer(data1, data2, 'Villar_2015')
    nbSigSnpInEnhancer(data1, data2, 'VISTA')
    nbSigSnpInEnhancer(data1, data2, 'FANTOM5')
    nbSigSnpInEnhancer(data1, data2, 'dbSUPER')
    nbSigSnpInEnhancer(resIntersection, data2, 'Villar_2015')
    nbSigSnpInEnhancer(resIntersection, data2, 'VISTA')
    nbSigSnpInEnhancer(resIntersection, data2, 'FANTOM5')
    nbSigSnpInEnhancer(resIntersection, data2, 'dbSUPER')

    g <- venn_diagram(data1, data2, pvalueSig)
    draw_plots(g, outfile2)
}

run()


