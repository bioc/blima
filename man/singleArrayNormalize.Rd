\name{singleArrayNormalize}
\alias{singleArrayNormalize}
\title{Bead level quantile normalization.}
\description{This function does quantile normalization of object beadLevelData from package beadarray.
Internal function not intended to direct use. Please use quantileNormalize.}
\usage{singleArrayNormalize(b, normalizationMod = NULL, channelNormalize = "Grn", 
    channelOutput = "qua", channelInclude = NULL, dst)}
\arguments{
  \item{b}{Object beadLevelData from package beadarray}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b.}
  \item{channelNormalize}{Name of channel to normalize.}
  \item{channelOutput}{Name of output normalized channel.}
  \item{channelInclude}{This field allows user to set channel with weights which have to be in {0,1}.
All zero weighted items are excluded from quantile normalization and the value asigned to such probes is a close to value which would be assigned to them if not being excluded. 
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).}
  \item{dst}{This field must be sorted. It is a distribution of values to assign to ports.
By default this distribution is computed using meanDistribution function.}
}



\author{Vojtěch Kulvait}





