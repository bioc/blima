\name{nonPositiveCorrectSingleArray}
\alias{nonPositiveCorrectSingleArray}
\title{Correct non positive}
\description{INTERNAL FUNCTION Correction for positive values only}
\usage{nonPositiveCorrectSingleArray(b, normalizationMod = NULL, channelCorrect = "GrnF", 
    channelBackgroundFilter = "bgf", channelAndVector = NULL)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelCorrect}{Name of channel to correct.}
  \item{channelBackgroundFilter}{Filtered beads will have weight 0 and non filtered weight 1.}
  \item{channelAndVector}{Represents vector to bitvise multiple to the channelBackgroundFilter vector.}
}



\author{Vojtech Kulvait}





