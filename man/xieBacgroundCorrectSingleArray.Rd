\name{xieBacgroundCorrectSingleArray}
\alias{xieBacgroundCorrectSingleArray}
\title{INTERNAL FUNCTION Xie background correct.}
\description{INTERNAL This function is not intended for direct use.
Background correction according to non parametric estimator in
Xie, Yang, Xinlei Wang, and Michael Story. 
"Statistical Methods of Background Correction for Illumina BeadArray Data."
Bioinformatics 25, no. 6 (March 15, 2009): 751-57. doi:10.1093/bioinformatics/btp040.
The method is applied on the bead level.}
\usage{xieBacgroundCorrectSingleArray(b, normalizationMod = NULL, negativeArrayAddresses, 
    channelCorrect, channelResult, channelInclude = NULL)}
\arguments{
  \item{b}{Single beadLevelData object.}
  \item{normalizationMod}{NULL for processing all spots in b. Otherwise specifies logical vector of the length equals to the number of arrays in b.}
  \item{negativeArrayAddresses}{Vector of addresses of negative control probes on array}
  \item{channelCorrect}{Slot to perform convolution correction.}
  \item{channelResult}{Result channel, if this channel exists it will be overwritten.}
  \item{channelInclude}{This field allows user to set channel with weights which have to be from {0,1}.
All zero weighted items are excluded from summarization.
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).}
}



\author{Vojtech Kulvait}





