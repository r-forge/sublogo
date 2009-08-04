\name{seqs.to.mat}
\Rdversion{1.1}
\alias{seqs.to.mat}
\title{seqs.to.mat}
\description{Calculate pairwise differences between sequences using a substitution matrix.}
\usage{
seqs.to.mat(seq.vec, subs.mat = NULL)
}
\arguments{
  \item{seq.vec}{DNA or protein sequences.}
  \item{subs.mat}{Substitution matrix with dimnames that match the letters used in the sequence data, or a character vector that specifies a common substitution matrix (as defined in the Biostrings package). NULL specifies that we will guess a suitable substitution matrix to match your input sequences (DNA=>identity, protein=>BLOSUM62).}
}
\details{
}
\value{The matrix of distances between each input sequence, with dimnames corresponding to either the sequences, or the sequence names (if they exist)}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{}