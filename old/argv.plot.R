## arguments: seqs,subsmat,title,subtitle,cutline,dendwidth%
## seqs - filename of sequences
## subsmat - standard name or filename of substitution matrix
## title - Text string title of the plot
## subtitle - Text string title of the plot
## cutline - numeric value of cutoff line for sublogos
## dendwidth - percent of width of middle dendrogram
a <- commandArgs(T)
print(a)
seq.file <- a[1]
subs.mat <- a[2]
tit <- a[3]
sub.tit <- a[4]
cutline <- as.numeric(a[5])
dend.width <- as.numeric(a[6])
plot.width <- as.numeric(a[7])
plot.height <- as.numeric(a[8])

## test for if named substitution matrix exists
data(list=subs.mat,package="Biostrings")
if(subs.mat%in%ls()){
  subs.mat <- get(subs.mat)
}else{
  ## so far only accepted format is result of R write.table
  subs.mat <- read.table(subs.mat)
}

## ghetto nonstandard seq input format
##seqs <- scan(seq.file,what="character",quiet=T)
##names(seqs) <- seqs

source("sublogo.dendrogram.R")

seqs <- read.fasta(seq.file)

outfile <- paste(seq.file,'pdf',sep='.')
pdf(outfile,w=plot.width,h=plot.height)
print(outfile)
print(plot.width)
print(plot.height)
print(seqs)
print(subs.mat)
m <- seqs.to.mat(seqs,subs.mat)
##print(m) too big for large datasets
sublogo.dendrogram(m,tit,sub.tit,seq.file,cutline,dend.width)
dev.off()
