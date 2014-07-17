# Functions for analysis of bead level data. 
# Author: Vojtech Kulvait
# Licence: GPL-3
##############################################################################

nonPositiveCorrect <- structure(function
        ###Correction for positive values only
        ##title<<Correct non positive
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        channelCorrect="GrnF",##<<Name of channel to correct.
        channelBackgroundFilter="bgf",##<<Filtered beads will have weight 0 and non filtered weight 1.
        channelAndVector=NULL##<<Represents vector to bitvise multiple to the channelBackgroundFilter vector.
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
    if(!is.null(channelAndVector))
    {
        channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelAndVector, "error")
    }
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, nonPositiveCorrectSingleArray(b[[i]], normalizationMod[[i]], channelCorrect, channelBackgroundFilter, channelAndVector));
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
        print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
    }
})

nonPositiveCorrectSingleArray <-function
        ###INTERNAL FUNCTION Correction for positive values only
        ##title<<Correct non positive
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        channelCorrect="GrnF",##<<Name of channel to correct.
        channelBackgroundFilter="bgf",##<<Filtered beads will have weight 0 and non filtered weight 1.
        channelAndVector=NULL##<<Represents vector to bitvise multiple to the channelBackgroundFilter vector.
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
        bgf = x[[i]][,channelCorrect]>0
        if(!is.null(channelAndVector))
        {
            bgf = bgf & x[[i]][,channelAndVector]
        }
        x = setWeights(x, wts=bgf, array=i, wtName=channelBackgroundFilter)
    }
    return(x)
}

bacgroundCorrect <- structure(function
        ###Background correction procedure selecting beads with background Intensity I_b |mean - I_b | > k*SD(I_bs) for exclusion.  
        ##title<<Data background correction.
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        channelBackground="GrnB",##<<Name of channel to normalize.
        k=3,##<<Parameter of method stringency (default is 3).
        channelBackgroundFilter="bgf",##<<Filtered beads will have weight 0 and non filtered weight 1. 
        channelAndVector=NULL##<<Represents vector to bitvise multiple to the channelBackgroundFilter vector.
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
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelBackground, "error")
    if(!is.null(channelAndVector))
    {
        channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelAndVector, "error")
    }
    
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, bacgroundCorrectSingleArray(b[[i]], normalizationMod[[i]], channelBackground, k, channelBackgroundFilter, channelAndVector));
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
        print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
    }
})

bacgroundCorrectSingleArray <-function
        ###Background correction procedure selecting beads with background Intensity I_b |mean - I_b | > k*SD(I_bs) for exclusion, internal.
        ##title<<Data background correction.
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        channelBackground="GrnB",##<<Name of channel to normalize.
        k=3,##<<Parameter of method stringency (default is 3).
        channelBackgroundFilter="bgf",##<<Filtered beads will have weight 0 and non filtered weight 1.
        channelAndVector=NULL##<<Represents vector to bitvise multiple to the channelBackgroundFilter vector.
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
        bgf = filterBg(x[[i]][,channelBackground], k)
        if(!is.null(channelAndVector))
        {
            bgf = bgf & x[[i]][,channelAndVector]
        }
        x = setWeights(x, wts=bgf, array=i, wtName=channelBackgroundFilter)
    }
    return(x)
}


filterBg <-function
        ###Background correction procedure selecting beads with background Intensity I_b |mean - I_b | > k*SD(I_bs) for exclusion, internal.  
        ##This function is not intended to direct use.
        ##title<<Bg correct vector
        (x,##<<Vector to correct
        k=3##<<Parameter of method stringency (default is 3).
)
{
    m = mean(x)
    std = sd(x)
    return(abs(x - m)<=k*std);
}

