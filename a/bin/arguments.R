res <- ''

s <- unlist(strsplit(res, ' '))
program  <- tail(unlist(strsplit(s[2], '/')), 1)
infile1  <- s[3]
infile2  <- s[4]
infile3  <- s[5]
infile4  <- s[6]
outfile1 <- s[7]

