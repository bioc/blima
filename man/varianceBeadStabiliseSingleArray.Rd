\name{varianceBeadStabiliseSingleArray}
\alias{varianceBeadStabiliseSingleArray}
\title{Bead level VST.}
\description{This function is not intended to direct use it takes single beadLevelData object and do bead level variance stabilisation.}
\usage{varianceBeadStabiliseSingleArray(b, normalizationMod = NULL, 
    quality = "qua", channelInclude = "bgf", channelOutput = "vst")}
\arguments{
  \item{b}{Object beadLevelData.}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b.}
  \item{quality}{Quality to analyze, default is "qua".}
  \item{channelInclude}{This field allows user to set channel with weights which have to be in {0,1}.
All zero weighted items are excluded from t-test. 
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to "bgf").}
  \item{channelOutput}{Output from VST.}
}



\author{Vojtech Kulvait}





