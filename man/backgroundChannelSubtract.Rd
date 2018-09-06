\name{backgroundChannelSubtract}
\alias{backgroundChannelSubtract}
\title{Background channel subtraction}
\description{Function to subtract one channel from another producing new channel. Standard graphic subtraction.}
\usage{backgroundChannelSubtract(b, normalizationMod = NULL, channelSubtractFrom = "GrnF", 
    channelSubtractWhat = "GrnB", channelResult = "Grn")}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for performing on all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelSubtractFrom}{Name of channel to subtract from.}
  \item{channelSubtractWhat}{Name of channel to subtract.}
  \item{channelResult}{Result channel, if this channel exists it will be overwritten.}
}



\author{Vojtěch Kulvait}




\examples{
if(require("blimaTestingData") && interactive())
{
    #To perform background correction on blimatesting object for two groups. Background correction is followed by correction for non positive data. Array spots out of selected groups will not be processed.
    data(blimatesting)
    #Prepare logical vectors corresponding to conditions A and E.
    groups1 = "A";
    groups2 = "E";
    sampleNames = list()
    c = list()
    for(i in 1:length(blimatesting))
    {
        p = pData(blimatesting[[i]]@experimentData$phenoData)
        c[[i]] = p$Group \%in\% c(groups1, groups2);
        sampleNames[[i]] = p$Name
    }
    #Background correction and quantile normalization followed by testing including log2TransformPositive transformation.
    blimatesting = bacgroundCorrect(blimatesting, normalizationMod=c, channelBackgroundFilter="bgf")
    blimatesting = nonPositiveCorrect(blimatesting, normalizationMod=c, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
}
}
