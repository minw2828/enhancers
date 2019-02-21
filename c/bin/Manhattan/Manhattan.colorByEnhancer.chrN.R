# /usr/bin/rscript

## Description:
## This script plots a manhattan plot over meta-analysis results.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 05 July 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf> <gender> <pheno> <infile> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
outfile1 <- args[3]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/metaAnalysis/snpName_pheno_lowerLimit_upperLimit_Z_p.14.2016-01-19.c.tab'
# infile2  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/snp/finalise.Villar_2015.snp.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/c/out/GSEA/Manhattan.colorByEnhancer/output.14.Villar_2015.2016-01-19.c.png'


# Get Support libraries and cleaning data #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile1, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile1, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

# Get global variables #
chrN <- as.numeric(f21(f22(outfile1, '/', 1), '[.]', 2))
db   <- f21(f22(outfile1, '/', 1), '[.]', 3)
VERYSIGNIFLINE <- 10e-20
GENOMEWIDELINE <- 10e-08
SUGGESTIVELINE <- 10e-05

# define functions #
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
    data <- fread(infile, header = FALSE)
    setnames(data, cns)
    return(data)
}

rawProcess1 <- function(data) {
    data[, chr := chrN]
    data[, pos := as.numeric(f21(snpName, ':', 2))]
    setnames(data, names(data)[length(names(data))], 'pos')
    return(data)
}

rawProcess2 <- function(data) {
    res <- data[chr == chrN]
    return(res)
}

combineDTs <- function(data1, data2) {
    setkey(data1, chr,pos)
    setkey(data2, chr,pos)
    mdt <- merge(data1, data2, all.x = TRUE)
    mdt[!is.na(source), boolenhancer := TRUE]
    mdt[is.na(source), boolenhancer := FALSE]
    return(mdt)
}

getColourList <- function() {
    dt <- data.table(colourName = c('verysignifline', 'genomewideline', 'suggestiveline'),
                     significance = c(VERYSIGNIFLINE, GENOMEWIDELINE, SUGGESTIVELINE))
    dt[, Colour := brewer.pal(12, "Paired")[as.factor(colourName)]]
    return(dt[])
}

plot_Manhattan <- function(data, colorList) {
    Title <- paste('Manhattan: Meta-analysis. Chromosome', chrN)
    SubTitle <- paste(colorList[, subTitle := paste(colourName, format(significance, scientific = TRUE), sep = ':')]$subTitle, collapse = ', ') 
    legend_name_fill <- "colour points"
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = pos, y = -log(pmeta), colour = boolenhancer, alpha = 0.8), size = 0.5) +
         geom_hline(mapping = aes(yintercept = -log(significance), colour = Colour), data = colorList, alpha = 0.7, linetype = 4, size = 0.5) +
         facet_wrap(~ pheno, nrow = 3) + 
         scale_x_discrete(breaks = as.numeric(levels(as.factor(data$chr))), labels = as.numeric(levels(as.factor(data$chr)))) +
         guides(colour = guide_legend(title = legend_name_fill, order = 1)) +
         scale_alpha(guide = "none") +
         ggtitle(bquote(atop(.(Title), atop(.(SubTitle))))) +
         labs(x = paste("chromosome", sprintf("%02d", chrN)), y = '-log(pvalue)') +
         theme(plot.title      = element_text(face = "bold", size = 80),
               axis.title.x    = element_text(size = 50),
               axis.title.y    = element_text(size = 50),
               strip.text.x    = element_text(size = 50),
               axis.text.x     = element_text(size = 50),
               axis.text.y     = element_text(size = 40),
               legend.text     = element_text(size = 40),
               legend.title    = element_text(size = 40))
    return(g)
}

zoom <- function(data, outfile) {
    return(data[which(snpPos >= ZOOMINSTART & snpPos <= ZOOMINEND)])
}

draw_plots1 <- function(outfile, g, Width = PLOTWIDTH, Height = PLOTHEIGHT) {
    Text <- paste(AIM, paste("-log", BASE, "(P)", sep = ""), paste(GENDER, PHENO, paste('Chr', sprintf("%02d", CHR), sep = ""), DB), sep = ": ")
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*0.5 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.5, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 50), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(g, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    dev.off()
}

# analysis #
run <- function() {
    data1 <- read_file(infile1)
    rawProcess1(data1)
    data2 <- fread(infile2, header = TRUE)
    dat2 <- rawProcess2(data2)
    data <- combineDTs(data1, dat2)
    colorList <- getColourList()
    g <- plot_Manhattan(data, colorList)
    f131(outfile1, g, 7200, 2880)
}

run()


