# /usr/bin/rscript

## Description:
## This module plots the properties of lesSNPs in each case.
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
## Rscript <module_name> <infile1> <infile2> <outfile1> <outfile12> <outfile13>


# Get input arguments #
args <- commandArgs(trailingOnly = TRUE)
infile1  <- args[1]
infile2  <- args[2]
infile3  <- args[3]
outfile1 <- args[4]
# infile1  <- '/tmp//1000bulls_v4_annotated_snps.tab'
# infile2  <- '/tmp//1000bulls_v4_annotated_indels.tab'
# infile3  <- '/tmp//CHROM_POS_db.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/categoryOfEnhancerVariants.1000BullGenomun4/piechart.functionalClass_enhancerVariants.2016-05-25.b.png'

# Get supports # 
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(unlist(strsplit(outfile1, '[.]')), 3)[1]
psf      <- tail(unlist(strsplit(outfile1, '[.]')), 2)[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

# analysis #
read_annoFile <- function(infile) {
    data <- fread(infile)
    setnames(data, '#CHROM', 'CHROM')
    data[CHROM == 'X', CHROM := 30]
    data[, CHROM := as.numeric(CHROM)]
    data[, POS := as.numeric(POS)]
    return(data)
}

read_enhancerFile <- function(infile) {
    cns <- unlist(strsplit(f21(f22((infile), '/', 1), '[.]', 1), '_'))
    data <- fread(infile)
    setnames(data, cns)
    return(unique(data))
}

mergeDataset <- function(data1, data2) {
    setkey(data1, CHROM,POS)
    setkey(data2, CHROM,POS)
    mdt <- merge(data1, data2, all.y = TRUE)
    return(mdt)
}

countFunctionalClassByDb <- function(data) {
    dt <- data.table(ddply(data, .(Functional_Class, db), nrow))
    setnames(dt, c(names(dt)[1:2], 'frequency'))
    dt[, Functional_Class := gsub('_', ' ', Functional_Class)]
    return(dt[])
}

piechart_countFunctionClass <- function(DT, Title) {
    bin_size <- 1
    Xlab  <- 'Count'
    legend_name_fill <- 'Functional Class'
    g <- ggplot(data = DT, aes(x = factor(1), y = frequency, fill = factor(Functional_Class))) +
         geom_bar(stat = "identity", width = bin_size) +
         ylim(0, 450000) + 
         coord_polar(theta = 'y') + 
#         facet_wrap(~db, ncol = 2) + 
         guides(fill = guide_legend(title = legend_name_fill)) +
         labs(title = Title, x = Xlab, y = '') +
         theme(plot.title = element_text(face = "bold", size = 20),
               axis.title.x = element_text(size = 15), 
               axis.title.y = element_text(size = 15),
               strip.text.x = element_text(size = 15), 
               legend.title = element_text(size = 15),
               legend.text  = element_text(size = 12),
               legend.position = "bottom") 
    return(g)
}

getDtPerAnnoType <- function(data1, data2, Title) {
    mdt <- mergeDataset(data1, data2)
    smdt <- mdt[, .(CHROM, POS, Functional_Class, Is_Fully_Known, db)]
    dt <- countFunctionalClassByDb(mdt)
    return(dt)
}

draw_plots12 <- function(Text, g1, g2, outfile, Width = 2400, Height = 1200) {
    # plot distbution ratio here is: Width : Height = 2400 : (960 * 4)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(1, 15), "null"))))
    grid.text(Text, gp = gpar(fontsize = 40), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    print(g1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(g2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
    dev.off()
}


## run ##
run <- function() {
    data1 <- read_annoFile(infile1)
    data2 <- read_annoFile(infile2)
    data3 <- read_enhancerFile(infile3)
    dt1 <- getDtPerAnnoType(data1, data3)
    dt2 <- getDtPerAnnoType(data2, data3)
    Text <- 'Count functional class of leading edge subset enhancer variants'
    g1 <- piechart_countFunctionClass(dt1[!is.na(Functional_Class)], '1000 Bull Genome Run4 SNP annotation')
    g2 <- piechart_countFunctionClass(dt2[!is.na(Functional_Class)], '1000 Bull Genome Run4 INDEL annotation')
    draw_plots12(Text, g1, g2, outfile1)
}

run()



