# /usr/bin/rscript

## Description:
## This script finds overlap between putative enhancers 
## (BLASTn search v.s. ChIP-Seq)
##
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 22 December 2015
##
## Date modified and reason:
##
## Execution:
## Rscript <module_name> <date> <psf> <infile> <outfile>


# get input arguments
args <- commandArgs(trailingOnly = TRUE)
date     <- args[1]
psf      <- args[2]
infile   <- args[3]
outfile  <- args[4] # venn diagram 
# date    <- '2016-05-25'
# psf     <- 'b' # project sub-folder
# infile  <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/snps/mergeAll.snps.tab'
# outfile <- '/group/dairy/Min/geno2pheno/analyses/2016-05-25/b/out/report/output.venn_diagram.snps.2016-05-25.b.png'

# define paths 
path_pre <- '/group/dairy/Min/geno2pheno'
path_me  <- file.path(path_pre, 'analyses', date, psf)
path_bin <- file.path(path_me, 'bin')
path_run <- file.path(path_me, 'run')
path_out <- file.path(path_me, 'out')
path_process <- file.path(path_out, 'process')
path_report  <- file.path(path_out, 'report')

# get support libraries and functions 
source(file.path(path_bin, 'libraries.R'))
source(file.path(path_bin, 'functions.R'))

# analysis
## read file ##
data <- fread(infile, header = TRUE)

data[, id := paste(chr, pos, sep = ":")]
dbs <- unique(data$source)

dt1 <- data[source == "VISTA"][][, "id", with=FALSE]
dt2 <- data[source == "FANTOM5"][, "id", with=FALSE]
dt3 <- data[source == "dbSUPER"][, "id", with=FALSE]
dt4 <- data[source == "Villar_2015"][, "id", with=FALSE]

fo12 <- Reduce(intersect, list(dt1$id, dt2$id))
fo13 <- Reduce(intersect, list(dt1$id, dt3$id))
fo14 <- Reduce(intersect, list(dt1$id, dt4$id))
fo23 <- Reduce(intersect, list(dt2$id, dt3$id))
fo24 <- Reduce(intersect, list(dt2$id, dt4$id))
fo34 <- Reduce(intersect, list(dt3$id, dt4$id))
fo123 <- Reduce(intersect, list(dt1$id, dt2$id, dt3$id))
fo124 <- Reduce(intersect, list(dt1$id, dt2$id, dt4$id))
fo134 <- Reduce(intersect, list(dt1$id, dt3$id, dt4$id))
fo234 <- Reduce(intersect, list(dt2$id, dt3$id, dt4$id))
fo1234 <- Reduce(intersect, list(dt1$id, dt2$id, dt3$id, dt4$id))

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
                            cex = rep(2, 15), 
                            cat.cex = rep(2, 4),
                            cat.col = colours)

draw_plots12 <- function(venn.plot, outfile, Width = 1200, Height = 960) {
    Text <- paste('Venn Diagram: Overlaps of putative bovine enhancer variants')
    # plot distbution ratio here is: Width : Height = 2400 : (960/5*1 + 960 + 960)
    Cairo(width = Width, height = Height, file = outfile, type = "png")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 5), "null"))))
    grid.text(Text, gp = gpar(fontsize = 30), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(grid.draw(venn.plot), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    dev.off()
}

draw_plots12(venn.plot, outfile)

