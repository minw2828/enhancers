# /usr/bin/rscript

## Description:
## This module plots up LD above FDR and colour dots by boolenhancer.
##
## Reference:
## 2005, Subramanian A. et al., Gene set enrichment analysis: A knowledge-based
## approach for interpreting genome-wide expression profiles, PNAS
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 26 April 2015
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <infile1> <infile2> <outfile> 


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1    <- args[1]
infile2    <- args[2]
outfile   <- args[3]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/bull/FY/14/Villar_2015/output.snpName_effect_pvalue_absEffectSorted_boolenhancer_pHitMiss_accSum.bull.FY.14.Villar_2015.2016-01-19.a.tab'
# infile2  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/bull/FY/14/Villar_2015/output.sigSnpA_sigSnpB_pearsonRcorrelation_pvalue.bull.FY.14.Villar_2015.0.0000001.2016-01-19.a.tab'
# outfile <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GSEA/bull/FY/14/Villar_2015/output.sigSnpA_sigSnpB_pearsonRcorrelation_pvalue.bull.FY.14.Villar_2015.0.0000001.0.01.2016-01-19.a.png'

# Get supports #
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

# Get global varaibles #
gender <- f21(f22(outfile, '/', 1), '[.]', 3)
pheno  <- f21(f22(outfile, '/', 1), '[.]', 4)
chrN   <- as.numeric(f21(f22(outfile, '/', 1), '[.]', 5))
db     <- f21(f22(outfile, '/', 1), '[.]', 6)
pThresholdGWAS <- as.numeric(paste(0, f21(f22(outfile, '/', 1), '[.]', 8), sep = '.'))
pearsonFDR     <- as.numeric(paste(0, f21(f22(outfile, '/', 1), '[.]', 10), sep = '.'))

# analysis #
## read file and give it colname ##
read_file <- function(infile) {
    cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 2), '_'))
    dat <- fread(infile, header = FALSE)
    setnames(dat, cns)
    return(dat)
}

filterFDR <- function(data2, pearsonFDR) {
    FDRthreshold <- compute.FDR(data2$pvalue, pearsonFDR)
    return(data2[which(pvalue <= FDRthreshold)])
}

rawProcess <- function(data1, data2) {
    sdat1 <- data1[, .(snpName, boolenhancer)]
    sdat2 <- filterFDR(data2, pearsonFDR)
    data <- get_boolenhancer(sdat1, sdat2)
    data[, ':=' (posA = as.numeric(f21(sigSnpA, ':', 2)), posB = as.numeric(f21(sigSnpB, ':', 2)))]
    return(data)
}

get_boolenhancer <- function(data1, data2) {
    setnames(data1, c('sigSnpA', 'boolenhancerA'))
    setkey(data1, sigSnpA)
    setkey(data2, sigSnpA)
    mdt <- merge(data1, data2)
    setnames(data1, c('sigSnpB', 'boolenhancerB'))
    setkey(data1, sigSnpB)
    setkey(mdt, sigSnpB)
    res <- merge(data1, mdt)
    return(res)
}

plot_LD_sigSNP <- function(data) {
    setkey(data, posA,posB)
    Title <- paste('LD between each SNPs in chromosome', chrN)
    Xlab <- 'snp position'
    Ylab <- 'snp position'
    g <- ggplot(data = data, mapping = aes(x = posA, y = posB, fill = pearsonRcorrelation)) +
         geom_tile(colour = "white", size = 0.1) +
         scale_fill_viridis(name = "pearson's r correlation", label = comma) +
         coord_equal() +
         labs(title = Title, x = Xlab, y = Ylab) +
         # element_text(family, face, colour, size)
         theme(plot.title = element_text(face = "bold", size = 25),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
               axis.text.x= element_text(size = 20),
               axis.text.y = element_text(size = 20))
    return(g)
}

draw_plots1 <- function(g, outfile, Width = 4800, Height = 2112) {
    Text <- paste(gender, pheno, chrN, db, sep = ', ')
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*0.5 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.5, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 30), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(g, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    dev.off()
}

## run ##
run <- function() {
    data1 <- read_file(infile1)
    data2 <- read_file(infile2)
    data <- rawProcess(data1, data2)
    g1 <- plot_LD_sigSNP(data)
    draw_plots1(g1, outfile)
}

run()



