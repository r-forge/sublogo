if(names(dev.cur())!="postscript"){ ## to avoid problems with R CMD check
sample.data <- c(
'QRTHLRD', 'QYVHLRS', 'QRAHLHS','QQSHLRA','QMAHLNS','QNTHLLQ',
'QNVHLTG','QRKDLRG','QNKDRAA','QPSNLHR','QKQNPIS','TSTHLQS',
'SWSGRRD')
library(sublogo)
sublogo(sample.data,cutline=30)

## simulation study: some seqs really are correlated, some are not,
## can we tell which ones are using these graphics?

## completely independent sequences
seqs <- replicate(40,paste(sample(dna.letters[2:5],13,rep=TRUE),collapse=''))
sublogo(seqs,cutline=10)
substr(seqs[1:15],1,1) <- "A"
substr(seqs[1:15],2,2) <- "G"
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
  lmat <- apply(m,1,function(p)sample(names(p),n,TRUE,p))
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

sublogo(seqs,cutline=6.7,main="Simulated sequences, with correlation between positions 1 and 10, but no correlation with position 4")

exsublogo <- function(S,subseq=NULL,...){
  seqs <- AlignedSeqs[[S]]
  if(!is.null(subseq))seqs <- substr(seqs,subseq[1],subseq[2])
  sublogo(seqs,...)
}
##debug(exsublogo)

data(AlignedSeqs)
exsublogo("zfp",cutline=30,main="Zinc finger protein recognition helix sequences, selected to bind triplet GGC")
exsublogo("cap.dna",dend.width=20,cutline=11.5,main="CAP promoters form a palindromic binding site motif",cex=0.75)
exsublogo("cap.protein",dend.width=25,cutline=75,cex=0.5,main="Helix-turn-helix motif from the Catabolite Activator Protein (CAP) transcription factor")
##exsublogo("globin",dend.width=20,cutline=90,cex=0.5,main="The end of the B helix through the beginning of the D helix of globins",subseq=c(61,81))
##exsublogo("prenyl",dend.width=30,cutline=160,main="Prenyltransferases (motif A)")
##exsublogo("splice",dend.width=20,cutline=27,main="Human splice sites on the intron/exon boundary")
}
