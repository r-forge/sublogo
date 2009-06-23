sample.data <- c(
'QRTHLRD', 'QYVHLRS', 'QRAHLHS','QQSHLRA','QMAHLNS','QNTHLLQ',
'QNVHLTG','QRKDLRG','QNKDRAA','QPSNLHR','QKQNPIS','TSTHLQS',
'SWSGRRD')
sample.data <- unique(sample.data)
pdf("sample.pdf",w=9.5,h=6.5)
source("sublogo.dendrogram.R")
sublogo(sample.data,cutline=30)
dev.off()

## simulation study: some seqs really are correlated, some are not,
## can we tell which ones are?

## completely independent sequences
seqs <- replicate(40,paste(sample(dna.letters[2:5],13,rep=T),collapse=''))
sublogo(seqs,cutline=10)
substr(seqs[1:15],1,1) <- "A"
substr(seqs[1:15],2,2) <- "G"
x11()
sublogo(seqs,cutline=10)

## better simulation: start out with a substitution matrix, set up a
## conditional probability link between positions i and j, simulate i
## according to substitution matrix, simulate j according to
## probability conditional on i, then add noise
