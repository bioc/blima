\name{backgroundChannelSubtractSingleArray}
\alias{backgroundChannelSubtractSingleArray}
\title{Background channel subtraction}
\description{INTERNAL FUNCTION Correction for positive values only}
\usage{backgroundChannelSubtractSingleArray(b, normalizationMod = NULL, 
    channelSubtractFrom = "GrnF", channelSubtractWhat = "GrnB", 
    channelResult = "Grn")}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelSubtractFrom}{Name of channel to subtract from.}
  \item{channelSubtractWhat}{Name of channel to subtract.}
  \item{channelResult}{Result channel, if this channel exists it will be overwritten.}
}



\author{Vojtech Kulvait}





