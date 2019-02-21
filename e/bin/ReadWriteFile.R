read_file <- function(infile) {
    info <- file.info(infile)
    if (is.na(info$isdir)) {
        return(NULL)
    } else if (info$size == 0) {
        return(NULL)
    } else {
        cns <- unlist(strsplit(f21(f22(infile, '/', 1), '[.]', 1), '_'))
        data <- fread(infile, header = FALSE)
        setnames(data, cns)
        return(data)
    }
}

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

deleteFiles <- function(files) {
    for (x in files) {
        if ( file.exists(x) ) {
            unlink(x)
        }
    }
}

writeLmSummary2File <- function(sf, outfile) {
    sink(outfile)
    print(sf)
    sink()
}


