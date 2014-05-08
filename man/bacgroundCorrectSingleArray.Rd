\name{bacgroundCorrectSingleArray}
\alias{bacgroundCorrectSingleArray}
\title{Data background correction.}
\description{Background correction procedure selecting beads with background Intensity I_b |mean - I_b | > k*SD(I_bs) for exclusion, internal.}
\usage{bacgroundCorrectSingleArray(b, normalizationMod = NULL, channelBackground = "GrnB", 
    k = 3, channelBackgroundFilter = "bgf", channelAndVector = NULL)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelBackground}{Name of channel to normalize.}
  \item{k}{Parameter of method stringency (default is 3).}
  \item{channelBackgroundFilter}{Filtered beads will have weight 0 and non filtered weight 1.}
  \item{channelAndVector}{Represents vector to bitvise multiple to the channelBackgroundFilter vector.}
}



\author{Vojtech Kulvait}





