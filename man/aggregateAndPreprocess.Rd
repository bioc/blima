\name{aggregateAndPreprocess}
\alias{aggregateAndPreprocess}
\title{Aggregate data}
\description{This function is not intended to direct use. It helps perform work of doProbeTTests function. For each probe it prints mean and sd of an quality.}
\usage{aggregateAndPreprocess(x, quality = "qua", transformation = NULL)}
\arguments{
  \item{x}{Two column matrix to agregate with columns "ProbeID" and quality.}
  \item{quality}{Quality to analyze, default is "qua".}
  \item{transformation}{Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TranformPositive.}
}



\author{Vojtěch Kulvait}





