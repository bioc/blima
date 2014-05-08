\name{filterBg}
\alias{filterBg}
\title{Bg correct vector}
\description{Background correction procedure selecting beads with background Intensity I_b |mean - I_b | > k*SD(I_bs) for exclusion, internal.  }
\usage{filterBg(x, k = 3)}
\arguments{
  \item{x}{Vector to correct}
  \item{k}{Parameter of method stringency (default is 3).}
}



\author{Vojtech Kulvait}





