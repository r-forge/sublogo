source("sublogo.dendrogram.R")
data(BLOSUM62,package="Biostrings")
sample.data <- c(
'QRTHLRD', 'QYVHLRS', 'QRAHLHS','QQSHLRA','QMAHLNS','QNTHLLQ',
'QNVHLTG','QRKDLRG','QNKDRAA','QPSNLHR','QKQNPIS','TSTHLQS',
'SWSGRRD')
sample.data <- unique(sample.data)
pdf("sample.pdf",w=9.5,h=6.5)
sublogo.dendrogram(seqs.to.mat(sample.data,BLOSUM62),
                 "Sample cluster",
                 "First test",
                 "./testseqs/zfp2",
                   cutline=30)
dev.off()
