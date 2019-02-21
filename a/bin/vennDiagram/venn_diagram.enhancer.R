# /usr/bin/rscript

## Description:
## This script finds overlap between putative enhancers 
## (BLASTn search v.s. ChIP-Seq)
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 11 March 2015
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name> <date> <psf> <infile> <outfile>


# get input arguments
args <- commandArgs(trailingOnly = TRUE)
date    <- args[1]
psf     <- args[2]
infile  <- args[3]
outfile <- args[4]
# date    <- '2016-05-25'
# psf     <- 'b' # project sub-folder
# infile  <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/enhancers/mergeAll.enhancers.csv'
# outfile <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/report/output.venn_diagram.enhancers.2016-05-25.b.png'

# define paths 
path_pre <- '/group/dairy/Min/geno2pheno'
path_me  <- file.path(path_pre, 'analyses', date, psf)
path_bin <- file.path(path_me, 'bin')
path_run <- file.path(path_me, 'run')
path_out <- file.path(path_me, 'out')

# get support libraries and functions 
source(file.path(path_bin, 'libraries.R'))
source(file.path(path_bin, 'functions.R'))

# analysis
## read file ##
data <- fread(infile)
cns <- c("chr", "start", "end", "source")
setnames(data, cns)

## convert data.table object into GenomicRanges object ##
dbs <- unique(data$source)
gr <- GRanges(seqnames = Rle(data$chr),
              ranges = IRanges(data$start, end = data$end),
              strand = Rle(rep("*", nrow(data))),
              source = data$source)

## example of code ## 
# overlap_venn <- function(dat1, dat2, name1, name2) {
#     fo <- findOverlaps(dat1, dat2)
#     p <- draw.pairwise.venn(area1      = queryLength(fo), 
#                             area2      = subjectLength(fo), 
#                             cross.area = length(fo), 
#                             category   = c(name1, name2),
#                             col        = rainbow(2),
#                             fill       = rainbow(2),
#                             alpha      = c(0.5, 0.5),
#                             label.col  = c(rep("black", 3)), )
#     return(p)
# }

## subsetByOverlaps(): extracts the elements in the query that overlap at
##                     least one element in the subject ##
## Reduce(): allow the GRanges objects to be collapsed across the whole of the
##           GRangesList object ##
fo12 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "FANTOM5"]))
fo13 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "dbSUPER"]))
fo14 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "Villar_2015"]))
fo23 <- Reduce(subsetByOverlaps, list(gr[gr$source == "FANTOM5"], gr[gr$source == "dbSUPER"]))
fo24 <- Reduce(subsetByOverlaps, list(gr[gr$source == "FANTOM5"], gr[gr$source == "Villar_2015"]))
fo34 <- Reduce(subsetByOverlaps, list(gr[gr$source == "dbSUPER"], gr[gr$source == "Villar_2015"]))
fo123 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "FANTOM5"], gr[gr$source == "dbSUPER"]))
fo124 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "FANTOM5"], gr[gr$source == "Villar_2015"]))
fo134 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "dbSUPER"], gr[gr$source == "Villar_2015"]))
fo234 <- Reduce(subsetByOverlaps, list(gr[gr$source == "FANTOM5"], gr[gr$source == "dbSUPER"], gr[gr$source == "Villar_2015"]))
fo1234 <- Reduce(subsetByOverlaps, list(gr[gr$source == "VISTA"], gr[gr$source == "FANTOM5"], gr[gr$source == "dbSUPER"], gr[gr$source == "Villar_2015"]))

colours <- c("maroon1", "gold1", "deepskyblue1", "darkseagreen1")
venn.plot <- draw.quad.venn(area1 = length(gr[gr$source == "VISTA"]), 
                            area2 = length(gr[gr$source == "FANTOM5"]), 
                            area3 = length(gr[gr$source == "dbSUPER"]), 
                            area4 = length(gr[gr$source == "Villar_2015"]), 
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
                            cex = rep(2, 15), 
                            cat.cex = rep(2, 4),
                            cat.col = colours)
 
draw_plots12 <- function(venn.plot, outfile, Width = 1200, Height = 960) {
    Text <- paste('Venn Diagram: Overlaps of putative bovine enhancers')
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 30), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(grid.draw(venn.plot), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    dev.off()
}

draw_plots12(venn.plot, outfile)

