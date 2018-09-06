\name{createSummarizedMatrix}
\alias{createSummarizedMatrix}
\title{Summarized value matrix.}
\description{This function creates summarized matrix of values of certain type.}
\usage{createSummarizedMatrix(b, spotsToProcess = NULL, quality = "qua", 
    channelInclude = "bgf", annotationTag = NULL)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{spotsToProcess}{NULL for processing all spots in b. Otherwise specifies logical vector of the length equals to the number of arrays in b.}
  \item{quality}{Quality to matrize.}
  \item{channelInclude}{This field allows user to set channel with weights which have to be from {0,1}.
All zero weighted items are excluded from summarization.
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to "bgf").}
  \item{annotationTag}{Tag from annotation file which to use in resulting matrix as colname.}
}



\author{Vojtěch Kulvait}




\examples{
if(require("blimaTestingData") && require("illuminaHumanv4.db") && interactive())
{
    #Create summarization of nonnormalized data from GrnF column.
    data(blimatesting)
    blimatesting = bacgroundCorrect(blimatesting, channelBackgroundFilter="bgf")
    blimatesting = nonPositiveCorrect(blimatesting, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    #Prepare logical vectors corresponding to conditions A(groups1Mod), E(groups2Mod) and both(processingMod).
    nonnormalized = createSummarizedMatrix(blimatesting, quality="GrnF", channelInclude="bgf",
            annotationTag="Name")
    head(nonnormalized)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
}
}
