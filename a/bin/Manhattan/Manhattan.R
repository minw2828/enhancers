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
infile1 <- args[1]
outfile <- args[2]
# infile1  <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/Manhattan_chrN_annotate/name_chr_pos_effect_pvalue.bull.FY.26.2016-01-19.a.tab'
# outfile1 <- '/group/dairy/Min/geno2pheno/analyses/2016-01-19/a/out/GWAS/Manhattan_chrN_annotate/output.Manhattan.-log10.bull.FY.14.2016-01-19.a.png'


# Get Support libraries and cleaning data
path_pre <- '/group/dairy/Min/geno2pheno'
date     <- tail(tail(unlist(strsplit(outfile, '[.]')), 3))[1]
psf      <- tail(tail(unlist(strsplit(outfile, '[.]')), 2))[1]
binpath  <- file.path(path_pre, 'analyses', date, psf, 'bin')
outpath  <- file.path(path_pre, 'analyses', date, psf, 'out')
source(file.path(binpath, 'functions.R'))
source(file.path(binpath, 'libraries.R'))
options(scipen = 999) # force not to use scienteific notation

# Get global variables #
base   <- as.numeric(f22(outfile, '[.]', 6)[1])
gender <- f22(outfile, '[.]', 5)[1]
pheno  <- f22(outfile, '[.]', 4)[1]
verysignifline <- 10e-20
genomewideline <- 10e-8
suggestiveline <- 10e-5

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

plot_Manhattan_annotate <- function(data) {
    SubTitle <- paste(paste('verysignifline', format(verysignifline, scientific = TRUE), sep = ':'), 
                      paste('genomewideline', format(genomewideline, scientific = TRUE), sep = ':'), 
                      paste('suggestiveline', format(suggestiveline, scientific = TRUE), sep = ':'), sep = '; ')
    Title <- paste('Manhattan', get_ylab(base), paste(gender, pheno, sep = ", "), sep = ": ")
    Xlab  <- "chromosomes"
    Ylab  <- get_ylab(base)
    legend_name_fill <- "chromosomes"
    g <- ggplot(data = data) +
         geom_point(mapping = aes(x = snpName, y = log_transform(pvalue, base), colour = as.factor(chr)), size = 1) +
         geom_hline(mapping = aes(yintercept = log_transform(verysignifline, base), alpha = 0.5), colour = "#FF0066", linetype = 4, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(genomewideline, base), alpha = 0.5), colour = "#00FF00", linetype = 2, size = 0.5) +
         geom_hline(mapping = aes(yintercept = log_transform(suggestiveline, base), alpha = 0.5), colour = "#0033FF", linetype = 1, size = 0.5) +
         scale_x_discrete(breaks = as.numeric(levels(as.factor(data$chr))), labels = as.numeric(levels(as.factor(data$chr)))) +
         guides(colour = guide_legend(title = legend_name_fill)) +
         #scale_colour_discrete(guide = "none") +
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
    return(g)
}

f131 <- function(file, plotg, Width = 4800, Height = 960) {
    Cairo(width = Width, height = Height, file = file, type = "png")
    print(plotg)
    dev.off()
}

appendLeadingZeros <- function(data, chr) {
    dat <- data[chr == chr]
    n <- max(nchar(dat$pos))
    dat[, snpName := paste(sprintf("%02d", chr), sprintf(paste("%0", n, "d", sep = ""), pos), sep = ":")]
    return(dat)
}

rawProcess <- function(infile) {
    data <- read_file(infile, unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_')))
    setkey(data, chr, pos)
    DATA <- rbindlist(lapply(seq(1, 30), function(x) appendLeadingZeros(data, x)))
    return(DATA)
} 

run <- function() {
    data <- rawProcess(infile1)
    g <- plot_Manhattan_annotate(data)
    f131(outfile, g)
}

run()

