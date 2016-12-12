\name{meanDistribution}
\alias{meanDistribution}
\title{Produce sorted double vector with mean distribution.}
\description{This function processes arrays in the object beadLevelData from package beadarray and returns sorted double vector.
The vector has length prvku. And the distribution of this vector is a "mean" of all distributions of distributionChannel quantity in arrays.
In case that probe numbers are different from prvku it does some averaging.}
\usage{meanDistribution(b, normalizationMod = NULL, distributionChannel = "Grn", 
    channelInclude = NULL, prvku)}
\arguments{
  \item{b}{Object beadLevelData from package beadarray or list of these objects}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical
vector of the length equals to the number of arrays in b or list of such vectors 
if b is a list of beadLevelData classes (defaults to NULL). }
  \item{distributionChannel}{Channel to do mean distribution from (defaults to "Grn").}
  \item{channelInclude}{This field allows user to set channel with weights which have to be in {0,1}.
All zero weighted items are excluded from quantile normalization and the value asigned to such probes is a close to value which would be assigned to them if not being excluded. 
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).}
  \item{prvku}{Number of items in a resulting double vector. Prvku must not be more than minimal number of indluded items in any distributionChannel.}
}



\author{Vojtěch Kulvait}





