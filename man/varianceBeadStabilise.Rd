\name{varianceBeadStabilise}
\alias{varianceBeadStabilise}
\title{Bead level VST.}
\description{This function does variance stabilising step on bead level.}
\usage{varianceBeadStabilise(b, normalizationMod = NULL, quality = "qua", 
    channelInclude = "bgf", channelOutput = "vst")}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for normalization of all input b. Otherwise specifies logical vector of the length equal to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes.}
  \item{quality}{Quality to analyze, default is "qua".}
  \item{channelInclude}{This field allows user to set channel with weights which have to be in {0,1}.
All zero weighted items are excluded from t-test.
You can turn this off by setting this NULL. This option is typically used together with background correction.}
  \item{channelOutput}{Output from VST.}
}



\author{Vojtech Kulvait}




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
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
}
}
