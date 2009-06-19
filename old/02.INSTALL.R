## Sublogo dendrograms
## Toby Dylan Hocking
## 25 February 2009

## To make sublogo dendrograms, you need to install:
## - R>=2.3 with CRAN packages grImport and gridBase.
install.packages(c("grImport","gridBase"))
## - Bioconductor package Biostrings.
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs="Biostrings")
## - Berkeley weblogo package (included, in weblogo/ subdirectory).
## - ImageMagick's convert command line image conversion program.

## Then for making sublogo dendrograms in plain R: cd to this directory,
## start the R interpreter, and try to follow the model in sample.R.

## For getting the web interface to work: you need to edit form.php and
## change the PATH= variable so that it includes the location of your R
## interpreter and convert programs. You need to have this directory
## accessible from a php-powered webserver, and then point your browser
## to form.php.

