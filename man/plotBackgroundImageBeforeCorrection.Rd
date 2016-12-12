\name{plotBackgroundImageBeforeCorrection}
\alias{plotBackgroundImageBeforeCorrection}
\title{Plot background image before correction}
\description{This function plots image of background distribution versus to foreground before background subtraction.}
\usage{plotBackgroundImageBeforeCorrection(b, index, channelForeground = "GrnF", 
    channelBackground = "GrnB", includePearson = FALSE)}
\arguments{
  \item{b}{Single beadLevelData object.}
  \item{index}{Index of spot to generate.}
  \item{channelForeground}{Name of channel of foreground.}
  \item{channelBackground}{Name of channel of background.}
  \item{includePearson}{Include Pearson corelation.}
}



\author{Vojtěch Kulvait}




\examples{
if(require("blimaTestingData") && interactive())
{
    #Write background images before correction. This function prints graph for condition D4. Call dev.off() to close.
    data(blimatesting)
    p = pData(blimatesting[[2]]@experimentData$phenoData)
    index = base::match("D4", p$Name)
    plotBackgroundImageBeforeCorrection(blimatesting[[2]], index)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
}
}
