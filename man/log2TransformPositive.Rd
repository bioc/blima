\name{log2TransformPositive}
\alias{log2TransformPositive}
\title{Log2 transform of numbers >1.}
\description{Transformation function are popular in beadarray package. Here this is similar concept. This function allow user to perform log transformation before doing t-tests.}
\usage{log2TransformPositive(x)}
\arguments{
  \item{x}{Number to transform.}
}

\value{This function returns logarithm of base 2 for numbers >=1 and zero for numbers <1.}

\author{Vojtech Kulvait}




\examples{
if(require("blimaTestingData") && require("illuminaHumanv4.db") && interactive())
{
    #To perform background correction, quantile normalization and then bead level t-test on log data run. Vst is not performed in this scheme. Top 10 probes is then printed according to certain measure.
    data(blimatesting)
    #Prepare logical vectors corresponding to conditions A(groups1Mod), E(groups2Mod) and both(c).
    groups1 = "A";
    groups2 = "E";
    sampleNames = list()
    groups1Mod = list()
    groups2Mod = list()
    c = list()
    for(i in 1:length(blimatesting))
    {
        p = pData(blimatesting[[i]]@experimentData$phenoData)
        groups1Mod[[i]] = p$Group \%in\% groups1;
        groups2Mod[[i]] = p$Group \%in\% groups2;
        c[[i]] = p$Group \%in\% c(groups1, groups2);
        sampleNames[[i]] = p$Name
    }
    #Background correction and quantile normalization followed by testing including log2TransformPositive transformation.
    blimatesting = bacgroundCorrect(blimatesting, normalizationMod =c, channelBackgroundFilter="bgf")
    blimatesting = nonPositiveCorrect(blimatesting, normalizationMod=c, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    blimatesting = quantileNormalize(blimatesting, normalizationMod=c, channelNormalize="GrnF", channelOutput="qua", channelInclude="bgf")
    beadTest <- doTTests(blimatesting, groups1Mod, groups2Mod,
            transformation=log2TransformPositive, quality="qua", channelInclude="bgf")
    symbol2address <- merge(toTable(illuminaHumanv4ARRAYADDRESS), toTable(illuminaHumanv4SYMBOLREANNOTATED))
    symbol2address <- symbol2address[,c("SymbolReannotated", "ArrayAddress") ]
    colnames(symbol2address) <- c("Symbol", "ArrayAddressID")
    beadTest = merge(beadTest, symbol2address, by.x="ProbeID", by.y="ArrayAddressID")
    beadTestID = beadTest[,c("ProbeID", "Symbol")]
    beadTestFC = abs(beadTest[,"mean1"]-beadTest[,"mean2"])
    beadTestP = beadTest[,"adjustedp"]
    beadTestMeasure = (1-beadTestP)*beadTestFC
    beadTest = cbind(beadTestID, beadTestMeasure)
    colnames(beadTest) <- c("ArrayAddressID", "Symbol", "difexBL")
    sortBL = sort(-beadTest[,"difexBL"], index.return=TRUE)$ix
    beadTop10 = beadTest[sortBL[1:10],]
    print(beadTop10)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData') and illuminaHumanv4.db by running biocLite('illuminaHumanv4.db').");
}
}
