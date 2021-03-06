# /usr/bin/rscript

## Description:
## This script plots a manhattan plot over GWAS results.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 07 March 2015
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
base         <- as.numeric(args[4])
zoomIn_start <- as.numeric(args[5])
zoomIn_end   <- as.numeric(args[6])
plot_width   <- as.numeric(args[7])
plot_height  <- as.numeric(args[8])
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/Manhattan_chrN_annotate/name_chr_pos_effect_pvalue.bull.FY.26.2016-01-19.a.tab'
# infile2  <- '/group/dairy/Min/geno2pheno/data/annotation/Raven2014BMC.table4.800K_snp.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/Manhattan_chrN_annotate/output.Manhattan.-log10.bull.FY.14.2016-01-19.a.png'
# base         <- as.numeric('10')
# zoomIn_start <- as.numeric('0')
# zoomIn_end   <- as.numeric('84648390')
# plot_width   <- as.numeric('2400')
# plot_height  <- as.numeric('960')

# Get Support libraries and cleaning data #
get_support <- function(outfile) {
    path_pre <- '/group/dairy/Min/geno2pheno'
    date     <- tail(unlist(strsplit(tail(unlist(strsplit(outfile, '/')), 1), '[.]')), 3)[1]
    psf      <- tail(unlist(strsplit(tail(unlist(strsplit(outfile, '/')), 1), '[.]')), 2)[1]
    binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
    source(file.path(binpath, 'functions.R'))
    source(file.path(binpath, 'libraries.R'))
    options(scipen = 999) # force not to use scienteific notation
    return(binpath)
}

# read file #
read_file <- function(infile, cns) {
    data <- fread(infile, header = FALSE)
    setnames(data, cns)
    return(data)
}

# transform pvalues #
log_transform <- function(pvalue, base = 10) {
    return(-log(pvalue, base))
}

# get plot ylab #
get_ylab <- function(base = 10) {
    if (base == 100) {
        return("-log100(p)")
    } else if (base == 10) {
        return("-log10(p)")
    } else if (base == 1) {
        return("-log(p)")
    }
}

# get plot title #
get_title <- function(outfile, base = 10) {
    aim    <- f21(f22(outfile, '/', 1), '[.]', 2)
    gender <- f22(f22(outfile, '/', 1), '[.]', 6)[1]
    pheno  <- f22(f22(outfile, '/', 1), '[.]', 5)[1]
    chrN   <- f22(f22(outfile, '/', 1), '[.]', 4)[1]
    Title <- paste(aim, get_ylab(base), paste(gender, pheno, chrN, sep = ", "), sep = ": ")
    return(Title)
}

# Create manhattan plot -log10 #
plot_Manhattan_annotate <- function(data, outfile, base, verysignifline, genomewideline, suggestiveline, plot_width, plot_height) {
    SubTitle <- paste(paste('verysignifline', format(verysignifline, scientific = TRUE), sep = ':'), 
                      paste('genomewideline', format(genomewideline, scientific = TRUE), sep = ':'), 
                      paste('suggestiveline', format(suggestiveline, scientific = TRUE), sep = ':'), sep = '; ')
    Title <- get_title(outfile, base)
    Xlab  <- "chromosome"
    Ylab  <- get_ylab(base)
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = snp_pos, y = log_transform(pvalue, base), colour = as.factor(snp_chr)), size = 1) +
         geom_hline(mapping = aes(yintercept = log_transform(verysignifline, base), alpha = 0.5), colour = "#FF0066", linetype = 4, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(genomewideline, base), alpha = 0.5), colour = "#00FF00", linetype = 2, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(suggestiveline, base), alpha = 0.5), colour = "#0033FF", linetype = 1, size = 0.5) +
         geom_text(mapping = aes(x = snp_pos, y = log_transform(pvalue, base), label = target), alpha = 0.7, colour = "#000000", size = 4, fontface = 2, hjust = 0, vjust = 0) +
         scale_x_discrete(breaks = as.numeric(levels(as.factor(data$snp_chr))), labels = as.numeric(levels(as.factor(data$snp_chr)))) +
         scale_colour_discrete(guide = "none") +
         scale_alpha(guide = "none") +
         ggtitle(bquote(atop(.(Title), atop(.(SubTitle))))) +
         labs(x = Xlab, y = Ylab) +
         theme(plot.title = element_text(face = "bold", size = 25),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
               axis.text.x = element_text(size = 15),
               axis.text.y = element_text(size = 15),
               legend.text = element_text(size = 15),
               legend.title = element_text(size = 15))
    f131(outfile, g, plot_width, plot_height)
}

zoom <- function(data, zoomIn_start = 0, zoomIn_end = Inf) {
    return(data[which(snp_pos >= zoomIn_start & snp_pos <= zoomIn_end)])
}

convert_phenoName <- function(phenoName) {
    if (tolower(phenoName) %in% c('fat', 'fat yield')) { 
        return('FY')
    } else if (tolower(phenoName) %in% c('milk', 'milk yield')) {
        return('MY')
    } else if (tolower(phenoName) %in% c('protein', 'milk yield')) {
        return('PY')
    } 
}

# analysis #
get_support(outfile1)

cns1  <- c("snp_names", "xaxis", "snp_chr", "snp_pos", "effect", "pvalue")
data1 <- read_file(infile1, cns1)
data1[, pheno := f22(f22(infile1, '/', 1), '[.]', 5)[1]]

data2 <- fread(infile2, header = TRUE)
data2[, snp_names := paste("Chr", chr, ":", pos, sep = "")]
data2[, pheno := unlist(lapply(trait, function(x) convert_phenoName(x)))]

data2[, target := gsub('NA', '', paste(within_gene, LD100kb, LD500kb, sep = ""))]
data1[, target := ""]

setkey(data1, snp_names,pheno)
setkey(data2, snp_names,pheno)
data1[data2, target := i.target]

verysignifline <- 10e-20
genomewideline <- 10e-8
suggestiveline <- 10e-5

plot_Manhattan_annotate(zoom(data1, zoomIn_start, zoomIn_end), outfile1, base, verysignifline, genomewideline, suggestiveline, plot_width, plot_height)


