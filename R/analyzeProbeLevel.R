# Functions for analysis of probe level data. 
# Author: Vojtech Kulvait
# Licence: GPL-3
###############################################################################

aggregateAndPreprocess <-function
        ###This function is not intended to direct use. It helps perform work of doProbeTTests function. For each probe it prints mean and sd of an quality.
        ##title<<Aggregate data
(x,##<<Two column matrix to agregate with columns "ProbeID" and quality.
quality="qua",##<<Quality to analyze, default is "qua".
transformation=NULL##<<Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TranformPositive.
)
{
    if(!is.null(transformation))
    {
        x[, quality] = sapply(x[, quality], transformation);
    }
    y=aggregate(x[, quality], list(x[,"ProbeID"] ), function(x) c(mean = mean(x), sd = sd(x)));
    y=cbind(y[,1], y[,2])
    colnames(y) <- c("ProbeID", "mean", "sd")
    return(y);
}

doProbeTTests <- structure(function
        ###This function does aggregated probe level t-tests on the data provided by the object beadLevelData from package beadarray.
        ##title<<T-test of bead level data.
        (b,##<<List of beadLevelData objects (or single object).
        c1,##<<List of logical vectors of data to assign to the first group (or single vector).
        c2,##<<List of logical vectors of data to assign to the second group (or single vector).
        quality="qua",##<<Quality to analyze, default is "qua".
        channelInclude="bgf",##<<This field allows user to set channel with weights which have to be  {0,1}. 
        ##All zero weighted items are excluded from t-test. 
        ##You can turn this off by setting this NULL. This option is typically used together with beadarray QC (defaults to "wts").
        correction="BY",##<<Multiple testing adjustment method as defined by p.adjust function, default is "BY".
        transformation=NULL##<<Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TranformPositive.
)
{
    waslist = checkIntegrity(b, "error")
    if(!waslist)
    {
        b = list(b)
        c1=list(c1);
        c2=list(c2);
    }
    checkIntegrityLogical(c1, b, "error");
    checkIntegrityLogical(c2, b, "error");
    channelExistsIntegrityWithLogicalVectorList(b, c1, quality, "error")
    channelExistsIntegrityWithLogicalVectorList(b, c2, quality, "error")
    if(!is.null(channelInclude))
    {
        channelExistsIntegrityWithLogicalVectorList(b, c1, channelInclude, "error")
        channelExistsIntegrityWithLogicalVectorList(b, c2, channelInclude, "error")
    }
    
    if(is.null(channelInclude))
    {
        nowts=TRUE;
    }else
    {
        nowts=FALSE;
    }
    a1 = matrix(0,0,2)
    a2 = matrix(0,0,2)
    for(i in 1:length(b))
    {
        for(j in 1:nrow(b[[i]]))
        {
            #print(c("[i,j]=",i, j,"."))
            if(nowts)
            {
                if(c1[[i]][[j]]==TRUE)
                {
                    a1 = rbind(a1, aggregateAndPreprocess(b[[i]][[j]][,c("ProbeID", quality)], quality, transformation)[,c("ProbeID","mean")])
                }
                if(c2[[i]][[j]]==TRUE)
                {
                    a2 = rbind(a2, aggregateAndPreprocess(b[[i]][[j]][,c("ProbeID", quality)], quality, transformation)[,c("ProbeID","mean")])
                }
            }else
            {
                if(c1[[i]][[j]]==TRUE)
                {
                    
                    a1 = rbind(a1, aggregateAndPreprocess(b[[i]][[j]][b[[i]][[j]][,channelInclude]==1,c("ProbeID", quality)], quality, transformation)[,c("ProbeID","mean")])
                }
                if(c2[[i]][[j]]==TRUE)
                {
                    a2 = rbind(a2, aggregateAndPreprocess(b[[i]][[j]][b[[i]][[j]][,channelInclude]==1,c("ProbeID", quality)], quality, transformation)[,c("ProbeID","mean")])
                }
            }
        }
    }
    s1=sort(a1[,"ProbeID"], index.return=TRUE)$ix
    s2=sort(a2[,"ProbeID"], index.return=TRUE)$ix
    a1 = a1[s1,]
    a2 = a2[s2,]
    probeIDs = unique(c(a1[,"ProbeID"],a2[,"ProbeID"]))
    results = matrix(0,length(probeIDs), 5)
    colnames(results) <- c("ProbeID", "p", "adjustedp", "mean1", "mean2")
    length1=length(a1[,"ProbeID"])
    length2=length(a2[,"ProbeID"])
    resultIndex = 1
    po1 = 1
    po2 = 1
    #print(c("Uz jsem tu! length1, length2", length1, length2))
    while(po1 <= length1 && po2 <= length2)
    {
        probe1 = a1[po1, "ProbeID"]
        probe2 = a2[po2, "ProbeID"]
        if(probe1 == probe2)
        {
            v1 = readToVector(a1, po1, length1, "mean")
            v2 = readToVector(a2, po2, length2, "mean")
            results[resultIndex, 1] = probe1;
            if(length(v1)==1 | length(v2)==1)
            {
                results[resultIndex, 2] = 1;
            }else
            {
                results[resultIndex, 2] = t.test(v1, v2, alternative="two.sided")$p.value
            }
            results[resultIndex, 4] = mean(v1)
            results[resultIndex, 5] = mean(v2)
            po1 = po1 + length(v1)
            po2 = po2 + length(v2)
        }else if(probe1 < probe2)
        {
            results[resultIndex, 1] = probe1;
            results[resultIndex, 2] = 100;
            po1 = getNextVector(a1, po1, length1);
        }else if(probe1 > probe2)
        {
            results[resultIndex, 1] = probe2;
            results[resultIndex, 2] = 100;
            po2 = getNextVector(a2, po2, length2);
        }
        resultIndex = resultIndex + 1;
        #print(c("[po1, po2]", po1, po2))
    }
    while(po1 <= length1)
    {
        probe1 = a1[po1,"ProbeID"]
        results[resultIndex, 1] = probe1;
        results[resultIndex, 2] = 100;
        po1 = getNextVector(a1, po1, length1);
        resultIndex = resultIndex + 1;
        #print(c("[po1]", po1))
    }
    while(po2 <= length2)
    {
        probe2 = a2[po2,"ProbeID"]
        results[resultIndex, 1] = probe2;
        results[resultIndex, 2] = 100;
        po2 = getNextVector(a2, po2, length2);
        resultIndex = resultIndex + 1;
        #print(c("[ po2]",  po2))
    }
    results[,3] = p.adjust(results[, 2], method=correction)
    return(results);
}, ex=function()
{
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
            groups1Mod[[i]] = p$Group %in% groups1;
            groups2Mod[[i]] = p$Group %in% groups2;
            processingMod[[i]] = p$Group %in% c(groups1, groups2);
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
})

