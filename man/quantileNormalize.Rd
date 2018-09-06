\name{quantileNormalize}
\alias{quantileNormalize}
\title{Bead level quantile normalization.}
\description{This function does quantile normalization of object beadLevelData from package beadarray.}
\usage{quantileNormalize(b, normalizationMod = NULL, channelNormalize = "Grn", 
    channelOutput = "qua", channelInclude = NULL, dst)}
\arguments{
  \item{b}{Object beadLevelData from package beadarray or list of these objects}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{channelNormalize}{Name of channel to normalize.}
  \item{channelOutput}{Name of output normalized channel.}
  \item{channelInclude}{This field allows user to set channel with weights which have to be in {0,1}.
All zero weighted items are excluded from quantile normalization and the value asigned to such probes is a close to value which would be assigned to them if not being excluded. 
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).}
  \item{dst}{User can specify sorted vector which represents distribution that should be assigned to items.}
}



\author{Vojtěch Kulvait}




\examples{
if(require("blimaTestingData") && interactive())
{
    #To perform background correction, variance stabilization and quantile normalization.
    data(blimatesting)
    #Prepare logical vectors corresponding to conditions A(groups1Mod), E(groups2Mod) and both(c).
    groups1 = "A";
    groups2 = "E";
    sampleNames = list()
    processingMod = list()
    for(i in 1:length(blimatesting))
    {
        p = pData(blimatesting[[i]]@experimentData$phenoData)
        processingMod[[i]] = p$Group \%in\% c(groups1, groups2);
        sampleNames[[i]] = p$Name
    }
    #Background correction and quantile normalization followed by testing including log2TransformPositive transformation.
    blimatesting = bacgroundCorrect(blimatesting, normalizationMod = processingMod, channelBackgroundFilter="bgf")
    blimatesting = nonPositiveCorrect(blimatesting, normalizationMod = processingMod, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    blimatesting = varianceBeadStabilise(blimatesting, normalizationMod = processingMod,
            quality="GrnF", channelInclude="bgf", channelOutput="vst")
    blimatesting = quantileNormalize(blimatesting, normalizationMod = processingMod,
            channelNormalize="vst", channelOutput="qua", channelInclude="bgf")
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
}
}
