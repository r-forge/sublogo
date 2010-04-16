\name{make.logo.ps}
\alias{make.logo.ps}
\title{make logo ps}
\description{Create a logo using weblogo, then read it in using grImport}
\usage{
make.logo.ps(helices, psbase)
}
\arguments{
  \item{helices}{Sequences to plot in the logo}
  \item{psbase}{Base filename for the logo postscript and xml files, should be the
full path name except for the trailing .ps}
}

\value{Grid picture grob as read using readPicture}

\author{Toby Dylan Hocking <toby.hocking@inria.fr>}





