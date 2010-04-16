\name{AlignedSeqs}
\Rdversion{1.1}
\alias{AlignedSeqs}
\docType{data}
\title{
Aligned sequence data ready for sublogo dendrogram plots.
}
\description{
This data set consists of a single named list: AlignedSeqs. Each element
is a character vector of aligned sequence data. For a given character
vector the sequence length of each element is constant. Thus you can use
these data sets with \code{sublogo} and \code{seqs.to.mat}.
}
\usage{data(AlignedSeqs)}
\format{
  The format is:
List of 6
 $ cap.dna    : Named chr [1:49] "attcgtgatagctgtcgtaaag" "ttttgttacctgcctctaactt" "aagtgtgacgccgtgcaaataa" "tgccgtgattatagacactttt" ...
  ..- attr(*, "names")= chr [1:49] "aldB-18->4" "ansB103->125" "araB1109->131" "araB2147->169" ...
 $ cap.protein: Named chr [1:101] "RGDIGNYLGLTVETISRLLG" "RGDIGNYLGLTVETISR---" "RGDIGNYLGLTVETISR---" "RGDIGNYLGLTVETISR---" ...
  ..- attr(*, "names")= chr [1:101] "RGDIGNYLGLTVETISRLLG" "RGDIGNYLGLTVETISR---" "RGDIGNYLGLTVETISR---" "RGDIGNYLGLTVETISR---" ...
 $ globin     : Named chr [1:56] "PIVDTGSVA-P------------------LSAAEKTKIRSAWAPVYSTY---ETSGVDILVKFFTSTPAAQEFFPKFKGLTT-----ADQLKKSA---DVRWHA-ERIINAVNDAVASMDDTEKMS-"| __truncated__ "PITDHGQPP-T------------------LSEGDKKAIRESWPQIYKNF---EQNSLAVLLEFLKKFPKAQDSFPKFSAKKS-------HLEQDP---AVKLQA-EVIINAVNHTIGLMDKEAAMK-"| __truncated__ "----------GS-----------------------DLVSGFWGKV--DA---HKIGGEALARLLVVYPWTQRYFTTFGNLGSADAIC-----HNA---KVLAHG-EKVLAAIGEGLKHPENLKAHY-"| __truncated__ "----------VN-----------------LTAAEKTQVTNLWGKV--NV---KELGGEALSRLLVVYPWTRRFFEHFGDLSTAEAVL-----HNA---KVLAHG-EKVLTSFGEGLKHLDNLKGTF-"| __truncated__ ...
  ..- attr(*, "names")= chr [1:56] "LampreyGLOBINV-SEALAMPREY" "HagfishGLOBINIII-ATLANTICHAGFISH" "FrogHEMOGLOBINBETACHAIN-EDIBLEFROG" "AfricanElephantHEMOGLOBINBETACHAIN-AFRICANELEPHANT" ...
 $ prenyl     : Named chr [1:14] "fktclevgev DTRGIYCALS iatllnILTE" "FGSCPHVdea HGGYTFCATA slailrsmdq" "GFCGRSNKlv DGCYSFWVGG saaileafgy" "DKPGAhsdfy HTNYCLLGLA vaessysctp" ...
  ..- attr(*, "names")= chr [1:14] "1DPR1$YEASTRAS-SPECIFICPRENYLTRANSFERASE(EC2.-.-.-)" "2DPR1$YEASTRAS-SPECIFICPRENYLTRANSFERASE(EC2.-.-.-)" "3DPR1$YEASTRAS-SPECIFICPRENYLTRANSFERASE(EC2.-.-.-)" "4DPR1$YEASTRAS-SPECIFICPRENYLTRANSFERASE(EC2.-.-.-)" ...
 $ splice     : Named chr [1:100] "ttctctgaaatatgaatttagACTGGTACTTATCATGGAG" "gcctgctttctcccctctcagGGACTTACAGTTTGAGATG" "cattgctgcttctttttttagGCATAAATTCTCGTGAACT" "aacttcctgtgtgttttgcagACAGCTGGATAGAAAACGA" ...
  ..- attr(*, "names")= chr [1:100] "19082_AF115399" "45328_AB000381" "45328_AB000381" "45330_AB001517" ...
 $ zfp        : Named chr [1:13] "DQGHRTR" "DVGHRSR" "ESGHLRR" "ESSKRKR" ...
  ..- attr(*, "names")= chr [1:13] "DQGHRTR" "DVGHRSR" "ESGHLRR" "ESSKRKR" ...
}
\source{
  These are the same data sets used in the original Weblogo paper:

  \url{http://weblogo.berkeley.edu}.
}
\examples{
## Can't be run since we use the mono font
## and the test initializes the postscript device without family="mono"
\dontrun{
data(AlignedSeqs)
library(sublogo)
sublogo(AlignedSeqs$zfp,
        cutline=30,
        main="Zinc finger protein recognition
          helix sequences, selected to bind triplet GGC")
}
}
\keyword{datasets}
