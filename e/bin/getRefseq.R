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

getTRefseq <- function(iRefseq, oRefseq) {
    # cannot convert to bostau8
    if (iRefseq == 'hg18') {
        if (oRefseq %in% c('hg19', 'hg38')) {
            TRefseq <- NA
        } else if (oRefseq %in% c('bostau6')) {
            TRefseq <- 'hg19'
        }
    } else if (iRefseq == 'mm9') {
        if (oRefseq %in% c('mm10', 'bostau6')) {
            TRefseq <- NA
        }
    } else if (iRefseq == 'canfam3') {
        if (oRefseq %in% c('bostau6')) {
            TRefseq <- 'canfam2'
        }
    } else if (iRefseq == 'mm10') {
        if (oRefseq %in% c('bostau6')) {
            TRefseq <- 'bostau7'
        }
    } else if (iRefseq == 'orycun2') {
        if (oRefseq %in% c('bostau6')) {
            TRefseq <- 'hg19'
       }
    } else if (iRefseq == 'rhemac2') {
        if (oRefseq %in% c('bostau6')) {
            TRefseq <- 'hg19'
        }
    } else {
        TRefseq <- NA
    }
    return(TRefseq)
}

