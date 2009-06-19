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
seqs.to.mat <- function(seq.vec,subs.mat){
  if(is.null(names(seq.vec)))names(seq.vec) <- seq.vec
  d <- toupper(gsub('[-.]',"*",seq.vec[unique(names(seq.vec))]))
  print(d)
  letters <- unique(c(unlist(strsplit(d,split='')),dna.letters))
  ##if dna alignment use simple identity matrix
  if(identical(sort(letters),sort(dna.letters))){
    subs.mat <- matrix(0,nrow=length(dna.letters),ncol=length(dna.letters),
                       dimnames=list(dna.letters,dna.letters))
    diag(subs.mat) <- 1
    subs.mat['*','*'] <- 0
  }
  print(subs.mat)
  N <- length(d)
  m <- matrix(0,nrow=N,ncol=N,dimnames=list(names(d),names(d)))
  for(i in 1:N)for(j in 1:i){
    seqs <- sapply(strsplit(c(d[i],d[j]),split=''),c)
    entry <- apply(seqs,1,function(x)subs.mat[x[1],x[2]])
    ## subscript out of bounds here usually means bad matrix
    m[i,j] <- m[j,i] <- -sum(entry)
  }
  m <- m-min(m)
  attr(m,'seqs') <- d
  m
}
##debug(seqs.to.mat)

## alternative vector format reading available with grImport + R>=2.3
make.logo.ps <- function(helices,psbase){
  epsfile <- paste(psbase,'eps',sep='.')
  psfile <- paste(psbase,'ps',sep='.')
  xmlfile <- paste(psfile,'xml',sep='.')
  ## to be replaced by berkeley weblogo eventually :
  #cmd <- paste('python helix_logo.py',"ATA",
  #             paste(helices,collapse=','),'>',psfile)
  seq.text <- paste(paste('>',helices,'\n',helices,sep=''),collapse='\n')
  write(seq.text,psbase)
  cmd <- paste("weblogo/seqlogo -c -F EPS -f",psbase,"|sed 's/^EndLine/%EndLine/'|sed 's/^EndLogo/%EndLogo/' >",psfile)#,'&& convert',epsfile,psfile)
  cat(cmd,'\n')
  system(cmd)
  PostScriptTrace(psfile,xmlfile)
  pic <- readPicture(xmlfile)
  #pic[-1:-41]
  pic
}
##debug(make.logo.ps)

subtitle <- function(st) mtext(st,line=-0.7,cex=1)
sublogo.dendrogram <- function(
  M,main,subtit,base,cutline=150,dend.width=30){
  hc <- hclust(as.dist(M),method="average")
  dend <- as.dendrogram(hc)
  fam <- cutree(hc,h=cutline)[labels(dend)] # order by plotting method
  nfams <- length(unique(fam))
  famtab <- data.frame(fam,y=1:length(fam))
  famtab$seq <- attr(M,'seqs')[rownames(famtab)]
  xrange <- c(0,1)
  yrange <- c(0,max(famtab[,'y'])+1)
  draw.box <- function(i){ # replace with logos
    subtab <- famtab[famtab[,'fam']==i,]
    grprange <- range(subtab[,'y'])
    rect(xrange[1],grprange[1]-0.25,xrange[2],grprange[2]+0.25)
  }
  draw.logo <- function(i){
    subrows <- famtab[,'fam']==i
    if(sum(subrows)>1){ # ignore family of size 1
      subtab <- famtab[subrows,]
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
      #addlogo(logo,xrange,grprange+c(-1,1)*0.25)
    }
  }

  bottomspace <- 0.6
  topspace <- 0.4
  ncex <- 1
  side.percents <- (100-dend.width)/2
  layout(matrix(1:3,ncol=3),c(side.percents,dend.width,side.percents))
  par(mai=c(bottomspace,0,topspace,0),cex=ncex)
  
  # Big summary logo on left
  plot(xrange,yrange,
       bty='n',xaxt='n',yaxt='n',ylab='',xlab='',type='n',yaxs='i')
  ##subtitle("Overall logo")
  biglogo <- make.logo.ps(famtab$seq,paste(base,'0',sep='.'))
  #addlogo(biglogo,xrange,range(famtab[,'y'])+c(0,1))
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
  
  # Dendrogram in middle
  ##par(mar=c(bottomspace,0,topspace,7.5))
  par(mai=c(bottomspace,0,topspace,
        max(strwidth(colnames(M),'inches',family='mono'))))
  par(family='mono')
  plot(dend,h=T,edgePar=list(lwd=2))
  par(family="")
  par(xpd=NA)
  segments(cutline,1,cutline,length(fam))
  
  # Title in the middle
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
  sapply(1:nfams,draw.logo)
  popViewport(3)
  
  hc
}
#debug(sublogo.dendrogram)


