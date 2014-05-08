\name{checkIntegrity}
\alias{checkIntegrity}
\title{Internal function}
\description{Check integrity of the list of beadLevelData objects or single beadLevelData object returns waslist.}
\usage{checkIntegrity(b, action = c("warn", "error"))}
\arguments{
  \item{b}{List of beadLevelData objects or single.}
  \item{action}{What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error.}
}

\value{Returns value if the object was list or not before calling this function.}

\author{Vojtech Kulvait}





