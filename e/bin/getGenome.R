getGenome <- function(data1, rowIndex, Refseq) {
    if (is.na(unique(data1[iRefseq == Refseq]$study))) {
        if (data1[iRefseq == Refseq]$iRefseq %in% c('bostau7')) {
             return(NULL)
        } else {
            Species <- data1[iRefseq == Refseq]$species
            refSeq  <- data1[iRefseq == Refseq]$refSeq
            ign <- paste('BSgenome', Species, 'UCSC', refSeq, 'masked', sep = '.')
            library(package = ign, character.only = TRUE)
            return(get(Species))
        }
    } else {
        if (data1[rowIndex][iRefseq == Refseq]$iRefseq %in% c('orycun2')) {
             return(NULL)
        } else {
            Species <- data1[rowIndex][iRefseq == Refseq]$species
            refSeq  <- data1[rowIndex][iRefseq == Refseq]$refSeq
            ign <- paste('BSgenome', Species, 'UCSC', refSeq, 'masked', sep = '.')
            library(package = ign, character.only = TRUE)
            return(get(Species))
        }
    }
}

