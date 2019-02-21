getSep <- function(infile) {
    Sep <- f22(infile, '[.]', 1)
    if (Sep %in% c('csv')) {
        return(',')
    } else if (Sep %in% c('bed', 'tab', 'sam')) {
        return('\t')
    } else {
        return(' ')
    }
}

read_file <- function(infile) {
    info <- file.info(infile)
    if (is.na(info$isdir)) {
        return(NULL)
    } else if (info$size == 0) {
        return(NULL)
    } else {
        Sep <- getSep(infile)
        cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
        data <- fread(infile, header = FALSE, sep = Sep)
        setnames(data, cns)
        return(data)
    }
}

write_file <- function(content, outfile, boolColname) {
    Sep <- getSep(outfile)
    write.table(content, outfile, append = FALSE, quote = FALSE, sep = Sep,
                row.names = FALSE, col.names = boolColname)
}

getGranges <- function(data) {
    if (is.null(data)) {
        return(NULL)
    } else {
        dat <- data[!(chr %in% c('chrUn'))]
        res <- as(dat, 'GRanges')
        seqlevelsStyle(res) <- 'UCSC'
        return(res)
    }
}

