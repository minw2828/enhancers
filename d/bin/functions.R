# /usr/bin/rscript

## Description:
## This script supports hist.R
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 24 November 2015
##
## Date modified and reason: 


# Analysis #
f11 <- function(file, plotg) {
    if (grep(".png", file)) {
        tmpname <<- sub(".png", ".pdf", x = file)
    } else if (grep(".jpeg", file)) {
        tmpname <<- sub(".jpeg", ".pdf", x = file)
    } else { warning("Unknown file format.") }
    ggsave(filename = tmpname, plot = plotg, width = 15, height = 12)
    ani.options('autobrowse' = FALSE)
    im.convert(files = tmpname, output = file, extra.opts="-density 150")
}

f12 <- function(content, path, qid, date, psf, fmt, sep, boo_colname) {
    outname <- paste('output', qid, date, psf, fmt, sep = ".")
    outfile <- file.path(path, outname)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = sep,
                row.names = FALSE, col.names = boo_colname)
}

f13 <- function(file, plotg) {
    Cairo(width = 1200, height = 960, file = file, type = "png")
    print(plotg)
    dev.off()
}

f131 <- function(file, plotg, Width = 1200, Height = 960) {
    Cairo(width = Width, height = Height, file = file, type = "png")
    print(plotg)
    dev.off()
}

f21 <- function(vec, sep, pos) {
    return(sapply(strsplit(vec, sep), function(x) x[pos]))
}

f22 <- function(vec, sep, pos) {
    return(sapply(strsplit(vec, sep), function(x) tail(x, pos)))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
   

