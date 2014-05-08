\name{vstFromLumi}
\alias{vstFromLumi}
\title{Function from LGPL lumi package 2.16.0 }
\description{This function is derived from copy and paste of lumi::vst function. Since lumi package has extensive imports I decided to hardcode this function to the blima instead of importing lumi package.}
\usage{vstFromLumi(u, std, nSupport = min(length(u), 500), backgroundStd = NULL, 
    lowCutoff = 1/3)}
\arguments{
  \item{u}{The mean of probe beads}
  \item{std}{The standard deviation of the probe beads}
  \item{nSupport}{Something for c3 guess.}
  \item{backgroundStd}{Estimate the background variance c3. Input should be variance according to article, not SD.}
  \item{lowCutoff}{Something for c3 guess.}
}


\references{\url{http://www.bioconductor.org/packages/release/bioc/html/lumi.html}}
\author{authors are Pan Du, Simon Lin, the function was edited by
Vojtech Kulvait}





