library(grImport)
library(gridBase)


## FASTA format
read.fasta <- function(infile){
  tmp <- sapply(strsplit(strsplit(paste('\n\n\n',paste(readLines(infile),collapse='\n'),sep=''),split='\n>')[[1]],'\n'),function(v)c(v[1],paste(v[2:length(v)],collapse='')))
  seqs <- tmp[2,]
  names(seqs) <- tmp[1,]
  blank <- names(seqs)==''
  names(seqs)[blank] <- seqs[blank]
  names(seqs) <- gsub(' ','',names(seqs))
  seqs[seqs!='']
}
##debug(read.fasta)

dna.letters <- c("*","A","T","G","C")
dna.identity <- matrix(0,nrow=length(dna.letters),ncol=length(dna.letters),
                       dimnames=list(dna.letters,dna.letters))
diag(dna.identity) <- 1
dna.identity['*','*'] <- 0
seqs.to.mat <- function(seq.vec,subs.mat=NULL){
  if(is.null(names(seq.vec)))names(seq.vec) <- seq.vec
  chars <- sapply(seq.vec,nchar)
  seqsum <- table(chars)
  if(length(seqsum)>1){
    print(as.data.frame(seqsum))
    i <- which.max(seqsum)
    cat("Sequences not of length ",names(seqsum)[i],":\n",sep="")
    rows <- as.integer(names(seqsum)[i])!=chars
    print(data.frame(chars,name=names(seq.vec),row.names=NULL)[rows,])
    stop("All input sequences must be of the same length.")
  }
  d <- toupper(gsub('[- .]',"*",seq.vec[unique(names(seq.vec))]))
  print(d)
  letters <- unique(c(unlist(strsplit(d,split='')),dna.letters))
  ##if dna alignment use simple identity matrix
  looks.like.dna <- identical(sort(letters),sort(dna.letters))
  if(is.null(subs.mat))
    subs.mat <- if(looks.like.dna)dna.identity else "BLOSUM62"
  if(mode(subs.mat)=="character"){
    ex <- substitute(data(M,package="Biostrings"),list(M="BLOSUM62"))
    eval(ex)
    subs.mat <- get(subs.mat)
  }
  print(subs.mat)
  N <- length(d)
  m <- matrix(0,nrow=N,ncol=N,dimnames=list(names(d),names(d)))
  for(i in 1:N)for(j in 1:i){
    seqs <- sapply(strsplit(c(d[i],d[j]),split=''),c)
    entry <- try(apply(seqs,1,function(x)subs.mat[x[1],x[2]]))
    if(class(entry)=="try-error"){
      print(seqs)
      stop("Sequence difference matrix construction failed.")
    }
    ## subscript out of bounds here usually means bad matrix
    m[i,j] <- m[j,i] <- -sum(entry)
  }
  m <- m-min(m)
  attr(m,'seqs') <- d
  m
}
##debug(seqs.to.mat)

## vector format reading available with grImport + R>=2.3
make.logo.ps <- function(helices,psbase){
  psfile <- paste(psbase,'ps',sep='.')
  xmlfile <- paste(psfile,'xml',sep='.')
  seq.text <- paste(paste('>',helices,'\n',helices,sep=''),collapse='\n')
  write(seq.text,psbase)
  cmd <- paste("PATH=/home/thocking/Desktop/sublogo/pkg/exec:$PATH seqlogo -c -F EPS -f",psbase,"|sed 's/^EndLine/%EndLine/'|sed 's/^EndLogo/%EndLogo/' >",psfile)
  cat(cmd,'\n')
  system(cmd)
  owd <- setwd(tempdir())
  PostScriptTrace(psfile,xmlfile)
  setwd(owd)
  pic <- readPicture(xmlfile)
  pic
}
##debug(make.logo.ps)

subtitle <- function(st) mtext(st,line=-0.7,cex=1)
sublogo.dendrogram <- function(
  M,main='',subtit=NULL,base=NULL,cutline=150,dend.width=30,cex=1){
  if(is.null(base))base <- tempfile()
  hc <- hclust(as.dist(M),method="average")
  dend <- as.dendrogram(hc)
  fam <- cutree(hc,h=cutline)[labels(dend)] # order by plotting method
  famids <- unique(fam)
  famtab <- data.frame(fam,y=1:length(fam))
  famtab$seq <- attr(M,'seqs')[rownames(famtab)]
  fam.nontriv <- sapply(famids,function(i)sum(famtab$fam==i))>1
  names(fam.nontriv) <- famids
  if(is.null(subtit))
    subtit <- paste(nrow(M),"sequences,",sum(fam.nontriv),"families")
  xrange <- c(0,1)
  yrange <- c(0,max(famtab[,'y'])+1)
  draw.box <- function(i){ # replace with logos
    subtab <- famtab[famtab[,'fam']==i,]
    grprange <- range(subtab[,'y'])
    rect(xrange[1],grprange[1]-0.25,xrange[2],grprange[2]+0.25)
  }
  draw.logo <- function(i){
    ## do not draw sublogo for a family of trivial size
    if(fam.nontriv[as.character(i)]){ 
      subtab <- famtab[famtab[,'fam']==i,]
      grprange <- range(subtab[,'y'])
      tmpfile <- paste(base,i,sep='.')
      logo <- make.logo.ps(subtab$seq,tmpfile)
      grid.picture(logo,
                   x=xrange[1],
                   y=unit(grprange[1]+0.25,'native'),
                   height=unit(diff(grprange)-0.5,'native'),
                   distort=T,
                   just=c(0,0),
                   fillText=T)
    }
  }

  bottomspace <- 0.6
  topspace <- 0.4
  ncex <- 1
  side.percents <- (100-dend.width)/2
  layout(matrix(1:3,ncol=3),c(side.percents,dend.width,side.percents))
  par(mai=c(bottomspace,0,topspace,0),cex=ncex)
  
  ## Big summary logo on left
  plot(xrange,yrange,
       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')
  ##subtitle("Overall logo")
  biglogo <- make.logo.ps(famtab$seq,paste(base,'0',sep='.'))
  vps <- baseViewports()
  pushViewport(vps$inner,vps$figure,vps$plot)
  grid.picture(biglogo,
               y=unit(0,'native'),
               height=unit(length(fam),'native'),
               just=c(0.5,0),
               exp=0,
               distort=T,
               fillText=T)
  popViewport(3)
  
  ## Dendrogram in middle
  par(mai=c(bottomspace,0,topspace,
        max(strwidth(colnames(M),'inches',
                     cex))/6*5),
      family='mono')
  par(xpd=NA)
  plot(dend,h=T,edgePar=list(lwd=2),nodePar=list(lab.cex=cex,pch=""))
  par(family="")
  segments(cutline,1,cutline,length(fam))
  ##axis(3,cutline,lty=0,line=0)
  
  ## Title in the middle
  title(main,line=0.5)
  ##subtitle("Dendrogram")
  subtitle(subtit)

  # Sublogos on the right side
  par(mai=c(bottomspace,0,topspace,0))
  plot(xrange,yrange,
       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')
  ##subtitle("Sublogos")
  vps <- baseViewports()
  pushViewport(vps$inner,vps$figure,vps$plot)
  sapply(famids,draw.logo)
  popViewport(3)
  
  hc
}
##debug(sublogo.dendrogram)

## shortcut function for common case of sequence data
sublogo <- function(seqs,mat=NULL,...)sublogo.dendrogram(seqs.to.mat(seqs,mat),...)
