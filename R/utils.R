# Help supporting functions. 
# Author: Vojtech Kulvait
# Licence: GPL-3
###############################################################################

doAction <- function
###Performs action of certain type
##title<<Internal function
(
message,##<<Text message.
action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error.
)
{
    action = match.arg(action)
    if(action == "returnText")
    {
        return(message);
    }else if(action == "warn")
    {
        warning(message);
        return(FALSE);
    }else if(action == "error")
    {
        stop(message);
    }
}

checkIntegrity <- function
###Check integrity of the list of beadLevelData objects or single beadLevelData object returns waslist.
##title<<Internal function
(
b,##<<List of beadLevelData objects or single.
action = c("warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(class(b)!="list")
    {
        b = list(b)
        waslist = FALSE;
    }else
    {
        waslist = TRUE;
    }
    checkIntegrityOfListOfBeadLevelDataObjects(b, action)
    return(waslist);
    ###Returns value if the object was list or not before calling this function.
}

checkIntegrityOfListOfBeadLevelDataObjects <- function
###Check integrity of the list of beadLevelData objects, internal.
##title<<Internal function
(listb,##<<List of beadLevelData objects.
action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(missing(listb))
    {
        return(doAction("Missing input arguments", action));
    }
    if(class(listb)=="list")
    {
        for(i in 1:length(listb))
        {
            res = checkIntegrityOfSingleBeadLevelDataObject(listb[[i]], "returnText")
            if(res != TRUE)
            {
                return(doAction(res, action));
            }
        }
    }else
    {
        return(doAction("This function should be used with list of objects of the type beadLevelData from package beadarray. Please ensure your object is compatible.", action))
    }
    return(TRUE);
}


checkIntegrityOfSingleBeadLevelDataObject <- function
        ###Check integrity of single beadLevelData object, internal.
        ##title<<Internal function
        (b,##<<beadLevelData object.
        action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(missing(b))
    {
        return(doAction("Missing input arguments", action));
    }
    if(!(class(b)=="beadLevelData"&&attributes(class(b))["package"]=="beadarray"))
    {
        return(doAction("This function should be used with object beadLevelData from package beadarray. Please ensure your object is compatible.", action))
    }
    return(TRUE);
}



checkIntegrityLogical <- function
###Check integrity of the list of logical objects, internal.
##title<<Internal function
(xx,##<<List of logical objects compatible with a list b.
b,##<<List of beadLevelData objects.
action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(missing(xx)||missing(b))
    {
        return(doAction("Missing input arguments", action));
    }
    if(class(b)=="list" && class(xx)=="list" && length(b) == length(xx))
    {
        for(i in 1:length(b))
        {
            if(!(class(xx[[i]])=="logical"&&length(xx[[i]])==nrow(b[[i]])))
            {
                return(doAction("Integrity error of logical vector.", action));
            }
        }
    }else
    {
        return(doAction("Integrity error of logical vector.", action));
    }
    return(TRUE);
}

singleCheckIntegrityLogicalVector <- function
        ###Check integrity of the logical object, internal.
        ##title<<Internal function
        (xx,##<<Logical object compatible with b.
        b,##<<Single beadLevelData object.
        action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(missing(xx)||missing(b))
    {
        return(doAction("Missing input arguments", action));
    }
    if(!(class(xx)=="logical"&&length(xx)==nrow(b)))
    {
        return(doAction("Integrity error of logical object.", action));
    }
    return(TRUE);
}

channelExistsIntegrityWithLogicalVectorList <- function
###Test existence of channel slot based on vector list
##title<<Internal function
(b,##<<List of beadLevelData objects.
spotsToCheck=NULL,##<<NULL for check all spots from b. Otherwise specifies logical vector of the length equals to the number of arrays in b with TRUE for checking.
slotToCheck,##<<Slot name to check
action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(missing(b)||missing(slotToCheck))
    {
        return(doAction("Missing input arguments", action));
    }
    if(!is.null(spotsToCheck))
    {
        res = checkIntegrityLogical(spotsToCheck, b, "returnText")
        if(res != TRUE)
        {
            return(doAction(res, action));
        }
    }
    for(i in 1:length(b))
    {
        res = singleChannelExistsIntegrityWithLogicalVector(b[[i]], spotsToCheck[[i]], slotToCheck, "returnText")
        if(res != TRUE)
        {
            return(doAction(res, action));
        }
    }    
    return(TRUE)
}

singleChannelExistsIntegrityWithLogicalVector <- function
        ###Test existence of channel slot based on logical list
        ##title<<Internal function
        (b,##<<single beadLevelData object
        spotsToCheck=NULL,##<<NULL for check all spots from b. Otherwise specifies logical vector of the length equals to the number of arrays in b with TRUE for checking.
        slotToCheck,##<<Slot name to check
        action = c("returnText", "warn", "error")##<<What type of action is required in case of invalid object structure. Either return text different from TRUE, warn or error. 
)
{
    action = match.arg(action)
    if(missing(b)||missing(slotToCheck))
    {
        return(doAction("Missing input arguments", action));
    }
    iteratorSet = 1:nrow(b);
    if(!is.null(spotsToCheck))
    {
        res = singleCheckIntegrityLogicalVector(spotsToCheck, b, "returnText")
        if(res != TRUE)
        {
            return(doAction(res, action));
        }
        iteratorSet = iteratorSet[spotsToCheck];
    }
    for(i in iteratorSet)
    {
        if(!(slotToCheck %in% colnames(b[[i]])))
        {
            return(doAction(sprintf("Slot %s is missing.", slotToCheck),action))
        }
    }
    return(TRUE);
}

log2TransformPositive <- structure(function
###Transformation function are popular in beadarray package. Here this is similar concept. This function allow user to perform log transformation before doing t-tests.
##title<<Log2 transform of numbers >1.
(x##<<Number to transform.
)
##value<<This function returns logarithm of base 2 for numbers >=1 and zero for numbers <1.
{
    if(x<1)
    {
        return(0);
    }else
    {
        return(log2(x))
    }
}, ex=function()
{
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
            groups1Mod[[i]] = p$Group %in% groups1;
            groups2Mod[[i]] = p$Group %in% groups2;
            c[[i]] = p$Group %in% c(groups1, groups2);
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
})
