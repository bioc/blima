\name{xieBacgroundCorrect}
\alias{xieBacgroundCorrect}
\title{Xie background correct.}
\description{Background correction according to non parametric estimator in
Xie, Yang, Xinlei Wang, and Michael Story. 
"Statistical Methods of Background Correction for Illumina BeadArray Data."
Bioinformatics 25, no. 6 (March 15, 2009): 751-57. doi:10.1093/bioinformatics/btp040.###The method is applied on the bead level.}
\usage{xieBacgroundCorrect(b, normalizationMod = NULL, negativeArrayAddresses, 
    channelCorrect, channelResult, channelInclude)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{normalizationMod}{NULL for processing all spots in b. Otherwise specifies logical vector of the length equals to the number of arrays in b.}
  \item{negativeArrayAddresses}{Vector of addresses of negative control probes on array}
  \item{channelCorrect}{Slot to perform convolution correction.}
  \item{channelResult}{Result channel, if this channel exists it will be overwritten.}
  \item{channelInclude}{This field allows user to set channel with weights which have to be from {0,1}.
All zero weighted items are excluded from summarization.
You can turn this off by setting this NULL. This option is typically used together with beadarray QC (defaults to "wts").}
}



\author{Vojtech Kulvait}




\examples{
if(require("blimaTestingData") && exists("annotationHumanHT12V4") && interactive())
{
    #Create vector of negative array addresses.
    negAdr = unique(annotationHumanHT12V4$Controls[annotationHumanHT12V4$Controls$Reporter_Group_Name=="negative", "Array_Address_Id"])
    #Create summarization of nonnormalized data from GrnF column.
    data(blimatesting)
    blimatesting = bacgroundCorrect(blimatesting, channelBackgroundFilter="bgf")
    blimatesting = nonPositiveCorrect(blimatesting, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    blimatesting = xieBacgroundCorrect(blimatesting, negativeArrayAddresses=negAdr, channelCorrect="GrnF", channelResult="GrnFXIE", channelInclude="bgf")
    #Prepare logical vectors corresponding to conditions A(groups1Mod), E(groups2Mod) and both(processingMod).
    xiecorrected = createSummarizedMatrix(blimatesting, quality="GrnFXIE", channelInclude="bgf",
            annotationTag="Name")
    head(xiecorrected)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData') and prepare annotationHumanHT12V4 object according to blimaTestingData manual.");
}
}
