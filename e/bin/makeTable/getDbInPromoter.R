getDbInPro <- function(pro, rgr3, nCov) {
    chr <- as.character(seqnames(pro))
    sgr3 <- rgr3[seqnames(rgr3) == chr]
    if (length(sgr3) != 0) {
        cr <- coverage(GRangesList(reduce(pro), sgr3))
        scr <- cr[[chr]]
        sscr <- slice(scr, lower = nCov)
    } else {
        sscr <- IRanges()
    }
    if (length(sscr) == 0) {
        res <- GRanges()
    } else {
        res <- GRanges(seqnames = chr, ranges(sscr))
    }
    return(res)
}

getRatioOfProInDb <- function(gr, gr3) {
    a <- sum(as.numeric(width(gr3)))
    b <- width(gr)
    res <- ifelse(is.na(a), 0, a/b)
    return(res)
}

getNbQtlInOutDb <- function(gr2, pro, his) {
    gr2.pro <- subsetByOverlaps(gr2, pro)
    gr2.his <- subsetByOverlaps(gr2, his) 
    nb2 <- length(gr2.pro)
    nb3 <- length(gr2.his)
    nb4 <- nb2 - nb3
    res <- data.table(nbQltInDb = nb3, nbQltOutDb = nb4)
    return(res)
}

