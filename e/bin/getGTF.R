getRefseq <- function(rf) {
    # hg18 == NCBI36 and mm9 == GRCm37 are both not available here
    if (rf == 'hg19') {
        refseq <- 'GRCh37'
    } else if (rf == 'mm10') {
        refseq <- 'GRCm38'
    } else if (rf == 'canfam3') {
        refseq <- 'CanFam3.1'
    } else if (rf == 'rhemac2') {
        refseq <- 'MMUL'
    } else if (rf == 'orycun2') {
        refseq <- 'OryCun2.0'
    } else if (rf == 'bostau6') {
        refseq <- 'UMD3.1'
    }
    return(refseq)
}

getGtfBySpeciesRefseq <- function(data1, ah, Refseq, tnms, cnms) {
    if (is.na(Refseq)) {
        return(NULL)
    } else if (!is.na(unique(data1[iRefseq == Refseq]$study)) & all(unique(data1[iRefseq == Refseq]$iRefseq) %in% c('hg18', 'mm9'))) {
        return(NULL)
    } else if (all(unique(data1[iRefseq == Refseq]$iRefseq) %in% c('canfam2', 'bostau7'))) {
        return(NULL)
    } else {
        RefSeq <- getRefseq(Refseq)
        Species <- unique(data1[iRefseq == Refseq]$SPECIES)
        ah.gr <- subset(ah, rdataclass == 'GRanges' & species == Species)
        gtfs <- query(ah.gr, c('gtf', RefSeq))
        n <- length(gtfs)
        gtf <- gtfs[[n]]
        seqlevelsStyle(gtf) <- 'UCSC'
        sgtf <- keepStandardChromosomes(gtf)
        res <- sgtf[sgtf$type %in% tnms, cnms]
        return(res)
    }
}

matchGtfTad <- function(gtf, gr) {
    if ((!is.null(gtf)) & (!is.null(gr))) {
        fo <- findOverlaps(gtf, gr)
        x <- gtf[queryHits(fo)]
        y <- gr[subjectHits(fo)]
        dt.x <- as.data.table(as.data.frame(x))
        dt.y <- as.data.table(as.data.frame(y))
        setnames(dt.x, paste(gsub('_', '', names(dt.x)), 'GTF', sep = ''))
        setnames(dt.y, paste(gsub('_', '', names(dt.y)), 'CTCF', sep = ''))
        res <- cbind(dt.x, dt.y)
        return(res)
    } else {
        return(NULL)
    }
}

