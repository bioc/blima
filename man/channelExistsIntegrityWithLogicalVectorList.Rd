\name{channelExistsIntegrityWithLogicalVectorList}
\alias{channelExistsIntegrityWithLogicalVectorList}
\title{Internal function}
\description{Test existence of channel slot based on vector list}
\usage{channelExistsIntegrityWithLogicalVectorList(b, spotsToCheck = NULL, 
    slotToCheck, action = c("returnText", "warn", "error"))}
\arguments{
  \item{b}{List of beadLevelData objects.}
  \item{spotsToCheck}{NULL for check all spots from b. Otherwise specifies logical vector of the length equals to the number of arrays in b with TRUE for checking.}
  \item{slotToCheck}{Slot name to check}
  \item{action}{What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error.}
}



\author{Vojtech Kulvait}





