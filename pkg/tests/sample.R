sample.data <- c(
'QRTHLRD', 'QYVHLRS', 'QRAHLHS','QQSHLRA','QMAHLNS','QNTHLLQ',
'QNVHLTG','QRKDLRG','QNKDRAA','QPSNLHR','QKQNPIS','TSTHLQS',
'SWSGRRD')
sample.data <- unique(sample.data)
pdf("sample.pdf",w=9.5,h=6.5)
source("sublogo.dendrogram.R")
sublogo(sample.data,cutline=30)
dev.off()
