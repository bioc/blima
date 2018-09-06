# Functions for analysis of probe level data. 
# Author: Vojtech Kulvait
# Licence: GPL-3
###############################################################################

xieBacgroundCorrect <- structure(function
###Background correction according to non parametric estimator in
###Xie, Yang, Xinlei Wang, and Michael Story. 
###"Statistical Methods of Background Correction for Illumina BeadArray Data."
###Bioinformatics 25, no. 6 (March 15, 2009): 751-57. doi:10.1093/bioinformatics/btp040.###The method is applied on the bead level.
##title<<Xie background correct.
(b,##<<List of beadLevelData objects (or single object).
normalizationMod=NULL,##<<NULL for processing all spots in b. Otherwise specifies logical vector of the length equals to the number of arrays in b. 
negativeArrayAddresses,##<<Vector of addresses of negative control probes on array
channelCorrect,##<<Slot to perform convolution correction.
channelResult,##<<Result channel, if this channel exists it will be overwritten.
channelInclude=NULL##<<This field allows user to set channel with weights which have to be from {0,1}.
##All zero weighted items are excluded from summarization.
##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).
)
{            
    waslist = checkIntegrity(b, "warn")
    if(!waslist)
    {
        b = list(b)
        if(!is.null(normalizationMod))
        {
            normalizationMod = list(normalizationMod)
        }
    }
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelCorrect, "error")
    if(!is.null(channelInclude))
    {
        channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelInclude, "error")
    }
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, xieBacgroundCorrectSingleArray(b[[i]], normalizationMod[[i]], negativeArrayAddresses, channelCorrect, channelResult, channelInclude));
    }
    if(waslist)
    {
        return(x);
    }else
    {
        return(x[[1]])
    }
}, ex=function()
{
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
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData') and prepare annotationHumanHT12V4 object according to blimaTestingData manual.");
    }
})
        
xieBacgroundCorrectSingleArray <- function
###INTERNAL This function is not intended for direct use.
###Background correction according to non parametric estimator in
###Xie, Yang, Xinlei Wang, and Michael Story. 
###"Statistical Methods of Background Correction for Illumina BeadArray Data."
###Bioinformatics 25, no. 6 (March 15, 2009): 751-57. doi:10.1093/bioinformatics/btp040.
###The method is applied on the bead level.
##title<<INTERNAL FUNCTION Xie background correct.
(b,##<<Single beadLevelData object.
normalizationMod=NULL,##<<NULL for processing all spots in b. Otherwise specifies logical vector of the length equals to the number of arrays in b. 
negativeArrayAddresses,##<<Vector of addresses of negative control probes on array
channelCorrect,##<<Slot to perform convolution correction.
channelResult,##<<Result channel, if this channel exists it will be overwritten.
channelInclude=NULL##<<This field allows user to set channel with weights which have to be from {0,1}.
##All zero weighted items are excluded from summarization.
##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn")
    x=b
    iteratorSet = 1:nrow(b);
    if(!is.null(normalizationMod))
    {
        iteratorSet = iteratorSet[normalizationMod];
    }
    for(i in iteratorSet)
    {
        toCorrectNeg = x[[i]][x[[i]][,channelInclude]==1 & x[[i]][,"ProbeID"] %in% negativeArrayAddresses,channelCorrect]
        toCorrectAll = x[[i]][x[[i]][,channelInclude]==1,channelCorrect]
        params <- nonParametricEstimator(toCorrectAll, toCorrectNeg);
        res = sapply(x[[i]][,channelCorrect], performXieCorrection, alpha=params$alpha, mu=params$mu, sigma=params$sigma)
        x = setWeights(x, wts=res, array=i, wtName=channelResult)
    }
    return(x)
}

nonParametricEstimator <- function
###INTERNAL This function is not intended for direct use.
###Background correction according to non parametric estimator in
###Xie, Yang, Xinlei Wang, and Michael Story. 
###"Statistical Methods of Background Correction for Illumina BeadArray Data."
###Bioinformatics 25, no. 6 (March 15, 2009): 751-57. doi:10.1093/bioinformatics/btp040.
###The method is applied on the bead level.
##title<<INTERNAL FUNCTION Xie background correct.
(toCorrectAll,#All convolution inputs.
toCorrectNeg#Normal distribution convolution inputs.
)
{
    params = list()
    meanAll = mean(toCorrectAll)
    meanNegative = mean(toCorrectNeg)
    params$alpha = meanAll - meanNegative
    params$mu = meanNegative
    params$sigma = sqrt(var(toCorrectNeg))
    return(params)
}

performXieCorrection <- function
###INTERNAL This function is not intended for direct use.
###Background correction according to non parametric estimator in
###Xie, Yang, Xinlei Wang, and Michael Story. 
###"Statistical Methods of Background Correction for Illumina BeadArray Data."
###Bioinformatics 25, no. 6 (March 15, 2009): 751-57. doi:10.1093/bioinformatics/btp040. ###The method is applied on the bead level.
##title<<INTERNAL FUNCTION Xie background correct.
(value, alpha, mu, sigma)
{
    a = value - mu - sigma^2/alpha
    val = a + sigma * dnorm(a/sigma) / pnorm(a/sigma);
    if(is.nan(val) || is.infinite(val))
    {
        val = 0
    }
    return(val)
}


backgroundChannelSubtract <- structure(function
                ###Function to subtract one channel from another producing new channel. Standard graphic subtraction.
                ##title<<Background channel subtraction
                (b,##<<List of beadLevelData objects (or single object).
                normalizationMod=NULL,##<<NULL for performing on all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
                channelSubtractFrom="GrnF",##<<Name of channel to subtract from.
                channelSubtractWhat="GrnB",##<<Name of channel to subtract.
                channelResult="Grn"##<<Result channel, if this channel exists it will be overwritten.
)
{
    waslist = checkIntegrity(b, "warn")
    if(!waslist)
    {
        b = list(b)
        if(!is.null(normalizationMod))
        {
            normalizationMod = list(normalizationMod)
        }
    }
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelSubtractFrom, "error")
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelSubtractWhat, "error")
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, backgroundChannelSubtractSingleArray(b[[i]], normalizationMod[[i]], channelSubtractFrom, channelSubtractWhat, channelResult));
    }
    if(waslist)
    {
        return(x);
    }else
    {
        return(x[[1]])
    }
}, ex=function()
{
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
            c[[i]] = p$Group %in% c(groups1, groups2);
            sampleNames[[i]] = p$Name
        }
        #Background correction and quantile normalization followed by testing including log2TransformPositive transformation.
        blimatesting = bacgroundCorrect(blimatesting, normalizationMod=c, channelBackgroundFilter="bgf")
        blimatesting = nonPositiveCorrect(blimatesting, normalizationMod=c, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    }else
    {
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})

backgroundChannelSubtractSingleArray <-function
        ###INTERNAL FUNCTION Correction for positive values only
        ##title<<Background channel subtraction
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        channelSubtractFrom="GrnF",##<<Name of channel to subtract from.
        channelSubtractWhat="GrnB",##<<Name of channel to subtract.
        channelResult="Grn"##<<Result channel, if this channel exists it will be overwritten.
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn")
    x=b
    iteratorSet = 1:nrow(b);
    if(!is.null(normalizationMod))
    {
        iteratorSet = iteratorSet[normalizationMod];
    }
    for(i in iteratorSet)
    {
        a = x[[i]][,channelSubtractFrom]-x[[i]][,channelSubtractWhat]
        x = setWeights(x, wts=a, array=i, wtName=channelResult)
    }
    return(x)
}


selectedChannelTransform <- structure(function
        ###Function to transform channel data.
        ##title<<Channel transformation
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for performing on all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        channelTransformFrom,##<<Name of channel to transform.
        channelResult,##<<Result channel, if this channel exists it will be overwritten.
        transformation=NULL##<<Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TranformPositive.
)
{
    waslist = checkIntegrity(b, "warn")
    if(is.null(transformation))
    {
        return(b);
    }
    if(!waslist)
    {
        b = list(b)
        if(!is.null(normalizationMod))
        {
            normalizationMod = list(normalizationMod)
        }
    }
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelTransformFrom, "error")
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, selectedChannelTransformSingleArray(b[[i]], normalizationMod[[i]], channelTransformFrom, channelResult, transformation));
    }
    if(waslist)
    {
        return(x);
    }else
    {
        return(x[[1]])
    }
}, ex=function()
{
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
            c[[i]] = p$Group %in% c(groups1, groups2);
            sampleNames[[i]] = p$Name
        }
        #Background correction and quantile normalization followed by testing including log2TransformPositive transformation.
        blimatesting = bacgroundCorrect(blimatesting, normalizationMod=c, channelBackgroundFilter="bgf")
        blimatesting = nonPositiveCorrect(blimatesting, normalizationMod=c, channelCorrect="GrnF",  channelBackgroundFilter="bgf", channelAndVector="bgf")
    }else
    {
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})

selectedChannelTransformSingleArray <- function
                ###Function to transform channel data.
                ##title<<Channel transformation
                (b,##<<List of beadLevelData objects (or single object).
                normalizationMod=NULL,##<<NULL for performing on all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
                channelTransformFrom,##<<Name of channel to transform.
                channelResult,##<<Result channel, if this channel exists it will be overwritten.
                transformation##<<Function of input data trasformation, default is NULL. Any function which for input value returns transformed value may be supplied. T-test then will be evaluated on transformed data, consider use log2TranformPositive.
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn")
    x=b
    iteratorSet = 1:nrow(b);
    if(!is.null(normalizationMod))
    {
        iteratorSet = iteratorSet[normalizationMod];
    }
    for(i in iteratorSet)
    {
        a = sapply(x[[i]][,channelTransformFrom], transformation)
        x = setWeights(x, wts=a, array=i, wtName=channelResult)
    }
    return(x)
}
