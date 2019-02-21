getValid <- function(tmp, genomeLen) {
    if (start(tmp) < 0)  {
        start(tmp) <- 1
    }
    if (end(tmp) > genomeLen) {
        end(tmp) <- genomeLen
    }
    return(tmp)
}

getPromoter <- function(gtf, pDis1, pDis2, genomeLen) {
    if(is.na(pDis1)& is.na(pDis2)) {
        tmp <- promoters(gtf)
    } else if (!is.na(pDis1) &  is.na(pDis2)) {
        tmp <- promoters(gtf, upstream = pDis1)
    } else if (is.na(pDis1) & !is.na(pDis2)) {
        tmp <- promoters(gtf, downstream = pDis2)
    } else if (!is.na(pDis1) & !is.na(pDis2)) {
        tmp <- promoters(gtf, upstream = pDis1, downstream = pDis2)
    } else {
        print('unmatched circumstances')
    } 
    res <- getValid(tmp, genomeLen)
    return(res)
}

