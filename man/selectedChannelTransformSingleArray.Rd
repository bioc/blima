\name{selectedChannelTransformSingleArray}
\alias{selectedChannelTransformSingleArray}
\title{Channel transformation}
\description{Function to transform channel data.}
\usage{selectedChannelTransformSingleArray(b, normalizationMod = NULL, 
    channelTransformFrom, channelResult, transformation)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for performing on all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelTransformFrom}{Name of channel to transform.}
  \item{channelResult}{Result channel, if this channel exists it will be overwritten.}
  \item{transformation}{Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TranformPositive.}
}



\author{Vojtěch Kulvait}





