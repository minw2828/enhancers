getData12 <- function(data1, data2) {
    data1[, `:=` (rowIndex = .I,
                  Factor   = paste(iRefseq, ct, sep = ':'))]
    setkey(data1, rowIndex)
    setkey(data2, rowIndex)
    mdt <- merge(data1, data2, all.y = TRUE)
    return(mdt)
}

getData34 <- function(data3, data4, pSigCTCF, ms) {
    dat3 <- data3[pvalue <= pSigCTCF & motifScore >= ms]
    data4[, Name := ifelse(Name == 30, 'X', as.character(Name))]
    setkey(dat3, RefSeq)
    setkey(data4, RefSeq)
    mdt <- merge(dat3, data4, allow.cartesian = TRUE)
    setnames(mdt, 'Name', 'chr')
    cns <- c('chr', setdiff(names(dat3), 'RefSeq'))
    res <- mdt[, cns, with = FALSE]
    return(res)
}

