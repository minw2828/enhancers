# /usr/bin/rscript

## Description:
## This script plots minor allele frequency by groups.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 15 March 2016
##
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf>


# Get input arguments
args <- commandArgs(trailingOnly = TRUE)
infile   <- args[1]
outfile  <- args[2]
# infile   <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/MAF/output.maf.jer.VISTA.2016-05-25.b.csv'
# outfile  <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/MAF/output.ProportionVariants_AlleleFrequency.jer.VISTA.2016-05-25.b.png'

# Get Support libraries and cleaning data
get_support <- function(outfile) {
    path_pre <- '/group/dairy/Min/geno2pheno'
    ofe      <- unlist(strsplit(unlist(strsplit(outfile, '/')), '[.]'))
    date     <- ofe[length(ofe)-2]
    psf      <- ofe[length(ofe)-1]
    binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
    source(file.path(binpath, 'functions.R'))
    source(file.path(binpath, 'libraries.R'))
    options(scipen = 999) # force not to use scienteific notation
    return(binpath)
}

# analysis
get_support(outfile)

# analysis
## read file ##
cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
data <- fread(infile, header = FALSE)
setnames(data, cns)

mt.data <- melt(data, id.vars = c("snpid", "boolenhancer"), measure.vars = c("freqAllele0", "freqAllele1"))

bin_size <- 0.001
nf <- nrow(mt.data[which(boolenhancer == FALSE)])
nt <- nrow(mt.data[which(boolenhancer == TRUE)])
f <- function(data, start, bin_size) {
    dat <- data[which(value >= start & value < start+bin_size),]
    dt <- as.data.table(ddply(dat, .(boolenhancer), nrow))
    cns <- c("boolenhancer", "count")
    setnames(dt, cns)
    dt[, proportion_of_variants := 0]
    dt[which(boolenhancer == FALSE)]$proportion_of_variants <- dt[which(boolenhancer == FALSE)]$count/nf
    dt[which(boolenhancer == TRUE)]$proportion_of_variants <- dt[which(boolenhancer == TRUE)]$count/nt
    dt[, allele_frequency_start := start]
    return(dt)
}

well <- lapply(seq(0, 1-bin_size, bin_size), function(x) f(mt.data, x, bin_size))
result <- rbindlist(well)

foo <- f22(f22(outfile, '/', 1), '[.]', 5)[1]
if (tolower(foo) == 'jer' | tolower(foo) == 'jersey') { 
    breed <- 'Jersey'
} else if (tolower(foo) == 'hol' | tolower(foo) == 'hos' | tolower(foo) == 'holstein') {
    breed <- 'Holstein'
} 
db    <- f22(f22(outfile, '/', 1), '[.]', 4)[1]
aim <- "Scatterplot: Proportion of Variants v.s. Allele Frequency"
Title <- paste(aim, paste("(Source: ", db, "; ", "Breed: ", breed, ")", sep = ""))
Xlab  <- paste("Allele Frequency;  bin size:", bin_size)
Ylab  <- "Proportion of Variants"
g <- ggplot(data = result) +
     geom_point(aes(x = allele_frequency_start, y = proportion_of_variants, colour = boolenhancer)) +
     facet_wrap(~ boolenhancer, ncol = 1) +
     labs(title = Title, x = Xlab, y = Ylab) + 
     theme(plot.title = element_text(face = "bold", size = 25),
           axis.title.x = element_text(size = 20),
           axis.title.y = element_text(size = 20),
           axis.text.x = element_text(size = 15),
           axis.text.y = element_text(size = 15),
           legend.text = element_text(size = 15),
           legend.title = element_text(size = 15))
f13(outfile, g)


