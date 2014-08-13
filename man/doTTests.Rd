\name{doTTests}
\alias{doTTests}
\title{T-test for bead (detector) level data.}
\description{This function does t-tests on the data provided by the object beadLevelData from package beadarray.}
\usage{doTTests(b, c1, c2, quality = "qua", channelInclude = "bgf", 
    correction = "BY", transformation = NULL)}
\arguments{
  \item{b}{List of beadLevelData objects (or single object).}
  \item{c1}{List of logical vectors of data to assign to the first group (or single vector).}
  \item{c2}{List of logical vectors of data to assign to the second group (or single vector).}
  \item{quality}{Quality to analyze, default is "qua".}
  \item{channelInclude}{This field allows user to set channel with weights which have to be  {0,1}.
All zero weighted items are excluded from t-test. 
You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to "bgf").}
  \item{correction}{Multiple testing adjustment method as defined by p.adjust function, default is "BY".}
  \item{transformation}{Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TransformPositive.}
}



\author{Vojtech Kulvait}




\examples{
if(require("blimaTestingData") && require("illuminaHumanv4.db") && interactive())
{
    #To perform background correction, variance stabilization and  quantile normalization then test on probe level, bead level and print top 10 results.
    data(blimatesting)
    #Prepare logical vectors corresponding to conditions A(groups1Mod), E(groups2Mod) and both(processingMod).
    groups1 = "A";
    groups2 = "E";
    sampleNames = list()
    groups1Mod = list()
    groups2Mod = list()
    processingMod = list()
    for(i in 1:length(blimatesting))
    {
        p = pData(blimatesting[[i]]@experimentData$phenoData)
        groups1Mod[[i]] = p$Group \%in\% groups1;
        groups2Mod[[i]] = p$Group \%in\% groups2;
        processingMod[[i]] = p$Group \%in\% c(groups1, groups2);
        sampleNames[[i]] = p$Name
    }
    #Background correction and quantile normalization followed by testing including log2TransformPositive transformation.
    blimatesting = bacgroundCorrect(blimatesting, normalizationMod =processingMod, channelBackgroundFilter="bgf")
    blimatesting = nonPositiveCorrect(blimatesting, normalizationMod=processingMod, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    blimatesting = varianceBeadStabilise(blimatesting, normalizationMod = processingMod, 
            quality="GrnF", channelInclude="bgf", channelOutput="vst")
    blimatesting = quantileNormalize(blimatesting, normalizationMod = processingMod, 
            channelNormalize="vst", channelOutput="qua", channelInclude="bgf")
    beadTest = doTTests(blimatesting, groups1Mod, groups2Mod, "qua", "bgf")
    probeTest = doProbeTTests(blimatesting, groups1Mod, groups2Mod, "qua", "bgf")
    adrToSymbol <- merge(toTable(illuminaHumanv4ARRAYADDRESS), toTable(illuminaHumanv4SYMBOLREANNOTATED))
    adrToSymbol <- adrToSymbol[,c("ArrayAddress", "SymbolReannotated") ]
    colnames(adrToSymbol) <- c("Array_Address_Id", "Symbol")
    probeTestID = probeTest[,"ProbeID"]
    beadTestID = beadTest[,"ProbeID"]
    probeTestFC = abs(probeTest[,"mean1"]-probeTest[,"mean2"])
    beadTestFC = abs(beadTest[,"mean1"]-beadTest[,"mean2"])
    probeTestP = probeTest[,"adjustedp"]
    beadTestP = beadTest[,"adjustedp"]
    probeTestMeasure = (1-probeTestP)*probeTestFC
    beadTestMeasure = (1-beadTestP)*beadTestFC
    probeTest = cbind(probeTestID, probeTestMeasure)
    beadTest = cbind(beadTestID, beadTestMeasure)
    colnames(probeTest) <- c("ArrayAddressID", "difexPL")
    colnames(beadTest) <- c("ArrayAddressID", "difexBL")
    tocmp <- merge(probeTest, beadTest)
    tocmp = merge(tocmp, adrToSymbol, by.x="ArrayAddressID", by.y="Array_Address_Id")
    tocmp = tocmp[, c("ArrayAddressID", "Symbol", "difexPL", "difexBL")]
    sortPL = sort(-tocmp[,"difexPL"], index.return=TRUE)$ix
    sortBL = sort(-tocmp[,"difexBL"], index.return=TRUE)$ix
    beadTop10 = tocmp[sortBL[1:10],]
    probeTop10 = tocmp[sortPL[1:10],]
    print(beadTop10)
    print(probeTop10)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData') and illuminaHumanv4.db by running biocLite('illuminaHumanv4.db').");
}
}
