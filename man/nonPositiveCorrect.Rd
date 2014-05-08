\name{nonPositiveCorrect}
\alias{nonPositiveCorrect}
\title{Correct non positive}
\description{Correction for positive values only}
\usage{nonPositiveCorrect(b, normalizationMod = NULL, channelCorrect = "GrnF", 
    channelBackgroundFilter = "bgf", channelAndVector = NULL)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelCorrect}{Name of channel to correct.}
  \item{channelBackgroundFilter}{Filtered beads will have weight 0 and non filtered weight 1.}
  \item{channelAndVector}{Represents vector to bitvise multiple to the channelBackgroundFilter vector.}
}



\author{Vojtech Kulvait}




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
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
}
}
