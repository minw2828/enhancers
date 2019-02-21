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


