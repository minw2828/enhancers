getSeqlens <- function(genome, Refseq) {
    if (!is.null(genome)) {
        tmp <- seqlengths(genome)
        seqlens <- tmp[names(tmp) %in% paste('chr', c(seq(1, 50), 'X', 'Y', 'M'), sep = '')]
    } else {
        if (Refseq == 'orycun2') {
            anno <- fread('/group/dairy/Min/geno2pheno/data/annotation/NCBI.refgen.infoTable.OryCun20.csv', sep = ',', header = TRUE)
            anno[, `:=` (chr    = paste('chr', Name, sep = ''), 
                         SizeBp = SizeMb * 1000000)]
            anno.sort <- anno[order(INSDC)]
        } else if (Refseq == 'bostau7') {
            anno <- fread('/group/dairy/Min/geno2pheno/data/annotation/NCBI.refgen.infoTable.bostau7.csv', sep = ',', header = TRUE)
            anno[, `:=` (chr    = paste('chr', Name, sep = ''),
                         SizeBp = Size)]
            anno.sort <- anno[order(INSDC)]
        }
        if (!('M' %in% anno.sort$Name)) {
            dt <- data.table(Loc = 'Nuc', Type = 'Chr', Name = 'M', RefSeq = NA, INSDC = NA, Size = NA, chr = 'chrM', SizeBp = NA)
            res <- rbind(anno.sort, dt)
            seqlens <- res$SizeBp
            names(seqlens) <- res$chr
        } else {
            seqlens <- anno.sort$SizeBp
            names(seqlens) <- anno.sort$chr
        }
    }
    return(seqlens)
}

getGaps <- function(genome, Refseq, gr) {
    anno <- getSeqlens(genome, Refseq)
    seqlengths(gr) <- anno[names(anno) %in% names(seqlengths(gr))]
    tmp  <- gaps(gr)
    res  <- tmp[strand(tmp) != '*']
    return(res)
}


