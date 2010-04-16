\name{sublogo.dendrogram}
\alias{sublogo.dendrogram}
\title{sublogo dendrogram}
\description{Main function for drawing sublogo dendrograms.}
\usage{
sublogo.dendrogram(M, main = "", subtit = NULL, base = tempfile(), cutline = 150, dend.width = 30, cex = 1)
}
\arguments{
  \item{M}{difference matrix as constructed by seqs.to.mat (although in
principle any object with a valid as.dist method could be used)}
  \item{main}{plot title}
  \item{subtit}{plot subtitle}
  \item{base}{base file name for temporary logo files}
  \item{cutline}{Distance for cutting the tree. Draw a sublogo for each
leaf. Normally you will plot once, then inspect the dendrogram to
determine which is a good value for cutline, then plot again using
your chosen cutline.}
  \item{dend.width}{Percent of the plot to be occupied by the dendrogram. The logos
will be displayed with equal widths on either side.}
  \item{cex}{character expansion factor for the dendrogram}
}

\value{The dendrogram from the call to hclust}

\author{Toby Dylan Hocking <toby.hocking@inria.fr>}





