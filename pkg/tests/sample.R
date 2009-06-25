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
good.rows <- function(x)rowSums(x>0.4)==2
rmat <- function(N=10){
  sm <- matrix(0)
  ## this generates an interesting starting substitution matrix, where
  ## there are at least 2 positions with freq>40%
  while(sum(gr <- good.rows(sm))<2){
    x <- matrix(rexp(N*4),N,4)
    sm <- x/rowSums(x)
  }
  print(cbind(gr,sm))
  colnames(sm) <- dna.letters[2:5]
  sm
}
rseq <- function(n,m,corr=0){
  lmat <- apply(m,1,function(p)sample(names(p),n,T,p))
  if(corr){
    good <- which(good.rows(m))
    p1 <- m[good[1],]
    g1 <- names(which(p1>0.4))
    p2 <- m[good[2],]
    g2 <- names(which(p2>0.4))
    favor <- function(l){
      i <- which(l==g1)
      if(length(i)==0)return(c(1,1,1,1))
      favored <- names(p2)==g2[i]
      (1-corr)/3*(!favored) + corr*favored
    }
    y <- t(sapply(lmat[,good[1]],favor))
    norm <- y/rowSums(y) ## easier to look at/compare
    lmat[,good[2]] <- apply(norm,1,function(v)sample(dna.letters[2:5],1,prob=v))    
  }
  apply(lmat,1,paste,collapse='')
}
##debug(rseq)
m <- rmat() ## make a substitution matrix
seqs <- rseq(50,m,0.95) ## randomly draw some sequences from that matrix

pdfname <- "../../poster/perfect-simulated-example.pdf"
pdf(pdfname,paper="a4r",h=0,w=0)
sublogo(seqs,cutline=6.7,main="Simulated sequences, with correlation between positions 1 and 10, but no correlation with position 4")
dev.off()
system(paste("xpdf",pdfname))

