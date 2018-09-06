\name{plotBackgroundImageAfterCorrection}
\alias{plotBackgroundImageAfterCorrection}
\title{Plot background image after correction}
\description{This function plots image of background distribution versus to foreground after background subtraction.}
\usage{plotBackgroundImageAfterCorrection(b, index, channelForeground = "GrnF", 
    channelBackground = "GrnB", SDMultiple = 3, includePearson = FALSE)}
\arguments{
  \item{b}{Single beadLevelData object.}
  \item{index}{Index of spot to generate.}
  \item{channelForeground}{Name of channel of foreground.}
  \item{channelBackground}{Name of channel of background.}
  \item{SDMultiple}{Correct on this level.}
  \item{includePearson}{Include Pearson corelation.}
}



\author{Vojtěch Kulvait}




\examples{
if(require("blimaTestingData") && interactive())
{
    #Write background images after correction. This function prints graph for condition D4. Call dev.off() to close.
    data(blimatesting)
    p = pData(blimatesting[[2]]@experimentData$phenoData)
    index = base::match("D4", p$Name)
    plotBackgroundImageAfterCorrection(blimatesting[[2]], index)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
}
}
