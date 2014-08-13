createSummarizedMatrix <- structure(function
        ###This function creates summarized matrix of values of certain type.
        ##title<<Summarized value matrix.
        (b,##<<List of beadLevelData objects (or single object).
        spotsToProcess=NULL,##<<NULL for processing all spots in b. Otherwise specifies logical vector of the length equals to the number of arrays in b. 
        quality="qua",##<<Quality to matrize.
        channelInclude="bgf",##<<This field allows user to set channel with weights which have to be from {0,1}.
        ##All zero weighted items are excluded from summarization.
        ##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to "bgf").
        annotationTag=NULL##<< Tag from annotation file which to use in resulting matrix as colname.
)
{
    waslist = checkIntegrity(b, "warn")
    if(!waslist)
    {
        b = list(b)
        if(!is.null(spotsToProcess))
        {
            spotsToProcess = list(spotsToProcess)
        }
    }
    checkIntegrityOfListOfBeadLevelDataObjects(b, "warn")
    channelExistsIntegrityWithLogicalVectorList(b, spotsToProcess, quality, "error")
    if(!is.null(channelInclude))
    {
        channelExistsIntegrityWithLogicalVectorList(b, spotsToProcess, channelInclude, "error")
    }
    output = NULL
    for(i in 1:length(b))
    {
        iteratorSet = 1:nrow(b[[i]])
        if(!is.null(spotsToProcess))
        {
            iteratorSet = iteratorSet[spotsToProcess[[i]]];
        }
        for(j in iteratorSet)
        {
            if(is.null(channelInclude))
            {
                setToPreprocess = b[[i]][[j]][,c("ProbeID", quality)]
            }else
            {
                setToPreprocess = b[[i]][[j]][b[[i]][[j]][,channelInclude]==1,c("ProbeID", quality)]
            }
            res = aggregateAndPreprocess(setToPreprocess, quality, NULL)[,c("ProbeID","mean")]
            if(is.null(annotationTag))
            {
                colnames(res)[2] = as.character(b[[i]]@sectionData$SampleGroup[j,])
            }else
            {
                pdata = pData(b[[i]]@experimentData$phenoData)
                colnames(res)[2] = pdata[j, annotationTag]
            }
            if(is.null(output))
            {
                output = res;
            }else
            {
                output = merge(output, res, by="ProbeID")
            }
        }
    }
    return(output)
}, ex=function()
{
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
        print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
    }
})
