# Functions for bead level quantile normalization of Illumina data processed by beadarray package. 
# Author: Vojtěch Kulvait
# Licence: GPL-3
###############################################################################

singleArrayNormalize <- function
###This function does quantile normalization of object beadLevelData from package beadarray.
###Internal function not intended to direct use. Please use quantileNormalize.
##title<<Bead level quantile normalization.
(b,##<<Object beadLevelData from package beadarray
normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b. 
channelNormalize="Grn",##<<Name of channel to normalize.
channelOutput="qua",##<<Name of output normalized channel.
channelInclude=NULL,##<<This field allows user to set channel with weights which have to be in {0,1}. 
##All zero weighted items are excluded from quantile normalization and the value asigned to such probes is a close to value which would be assigned to them if not being excluded. 
##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).
dst##<<This field must be sorted. It is a distribution of values to assign to ports.
##By default this distribution is computed using meanDistribution function.
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn")
    if(missing(dst))
        stop("Field dst is mandatory");
    x=b
    if(is.null(channelInclude))
    {
        nowts=TRUE;
    }else
    {
        nowts=FALSE;
    }
    
    bdc = integer();
    sw = integer();
    iteratorSet = 1:nrow(b);
    if(!is.null(normalizationMod))
    {
        iteratorSet = iteratorSet[normalizationMod];
    }
    for(i in iteratorSet)
    {
        bdc[i] = nrow(b[[i]])
        if(nowts)
        {
            sw[i] = bdc[i];
        }else
        {
            sw[i] = sum(b[[i]][,channelInclude])#\in {0,1}
        }
    }
    prvku = length(dst);
    
    for(i in iteratorSet)
    {
        q = double(bdc[i])
        sorted = sort(b[[i]][,channelNormalize], index.return=TRUE)$ix
        cnt = 0;
        if(nowts)
        {
            for(j in 1:bdc[i])
            {
                cnt = cnt + 1;
                pos = round((sw[i]-prvku+(prvku-1)*cnt)/(sw[i]-1))
                q[sorted[j]] = dst[pos]
            }
        }else
        {
            sortedweights = b[[i]][, channelInclude][sorted]
            for(j in 1:bdc[i])
            {
                if(sortedweights[j]==1)
                {
                    cnt = cnt + 1;
                    pos = round((sw[i]-prvku+(prvku-1)*cnt)/(sw[i]-1))
                    q[sorted[j]] = dst[pos]
                }else
                {
                    if(cnt==0)
                    {
                        q[sorted[j]] = dst[1]
                    }else
                    {
                        pos = round((sw[i]-prvku+(prvku-1)*cnt)/(sw[i]-1))
                        q[sorted[j]] = dst[pos]
                    }
                }
            }
        }
        if(cnt!=sw[i])
            stop("ERROR HERE, PLEASE REPORT BUG!")
        x = setWeights(x, wts=q, array=i, wtName=channelOutput)
    }
    return(x)
}

singleNumberOfDistributionElements <-function
###Internal function
##title<<Internal
(b,##<<Object beadLevelData from package beadarray
normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
channelInclude=NULL##<This field allows user to set channel with weights which have to be in {0,1}. 
)
{
    if(is.null(channelInclude))
    {
        nowts=TRUE;
    }else
    {
        nowts=FALSE;
    }
    #object inicialization
    sumWeights <- integer();
    
    iteratorSet = 1:nrow(b);
    if(!is.null(normalizationMod))
    {
        iteratorSet = iteratorSet[normalizationMod];
    }
    for(i in iteratorSet)
    {
        if(nowts)
        {
            sumWeights = c(sumWeights, nrow(b[[i]]));
        }else
        {
            sumWeights = c(sumWeights, sum(b[[i]][,channelInclude]))#\in {0,1}
        }
    }
    return(min(sumWeights))
}

numberOfDistributionElements <-function
###Internal function
##title<<Internal
(b,##<<Object beadLevelData from package beadarray or list of these objects
normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
channelInclude=NULL##<This field allows user to set channel with weights which have to be in {0,1}. 
)
{
    if(class(b)=="list")
    {
        sumWeights <- integer();
        if(is.null(normalizationMod))
        {
            for(i in 1:length(b))
            {
                sumWeights[i] = singleNumberOfDistributionElements(b[[i]], NULL, channelInclude)
            }
        }else
        {
            for(i in 1:length(b))
            {
                sumWeights[i] = singleNumberOfDistributionElements(b[[i]], normalizationMod[[i]], channelInclude)
            }
        }
        return(min(sumWeights));
    }else
    {
        return(singleNumberOfDistributionElements(b, normalizationMod, channelInclude));
    }
}

quantileNormalize <- structure(function
###This function does quantile normalization of object beadLevelData from package beadarray.
##title<<Bead level quantile normalization.
(b,##<<Object beadLevelData from package beadarray or list of these objects
normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
channelNormalize="Grn",##<<Name of channel to normalize.
channelOutput="qua",##<<Name of output normalized channel.
channelInclude=NULL,##<<This field allows user to set channel with weights which have to be in {0,1}. 
##All zero weighted items are excluded from quantile normalization and the value asigned to such probes is a close to value which would be assigned to them if not being excluded. 
##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).
dst##<<User can specify sorted vector which represents distribution that should be assigned to items.
#By default this distribution is computed using meanDistribution function.
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
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelNormalize, "error")
    if(!is.null(channelInclude))
    {
        channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelInclude, "error")
    }

    if(missing(dst) || is.null(dst))
    {
        prvku <- numberOfDistributionElements(b, normalizationMod, channelInclude);
        dst <- meanDistribution(b, normalizationMod, channelNormalize, channelInclude, prvku);
    }
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, singleArrayNormalize(b[[i]],normalizationMod[[i]],channelNormalize,channelOutput,channelInclude,dst));
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
            processingMod[[i]] = p$Group %in% c(groups1, groups2);
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
})

initMeanDistribution <-function
###This is internal function not intended to direct use which initializes mean distribution.
(srt,##<<vector of sorted values
prvku ##<<number of items in meanDistribution
)
{
    if(prvku>length(srt))
    {
        stop("Length of meanDistribution is less or equal length of particullar array!");
    }
#    meanDistribution <- double(prvku);
#    vsz = length(srt)
#    prev = 1;
#    prvek = srt[1];
#    if(vsz>1)
#    {
#        for(j in 2:vsz)
#        {
#            pos = round((vsz-prvku+(prvku-1)*j)/(vsz-1))
#            if(pos==prev)
#            {
#                prvek=c(prvek, srt[j])
#            }else
#            {
#                meanDistribution[pos-1] = mean(prvek);
#                prvek = srt[j];
#                prev = pos;
#            }
#        }
#    }
#    meanDistribution[prev] = mean(prvek);
#    return(meanDistribution);
	return(interpolateSortedVector(srt, prvku))
	
}

updateMeanDistribution<-function
###This is internal function not intended to direct use. Updates mean distribution.
(meanDistribution, #semicreated meanDistribution
srt, ##<<vector of sorted values
arraysUsed ##<<number of arrays allready used to create distribution
)
{
    if(length(meanDistribution)>length(srt))
    {
        stop("Length of meanDistribution is less or equal length of particullar array!");
    }
#    if(length(meanDistribution)<3)
#    {
#        warning("DEBUG:Something weird with meanDistribution!")
#    }
#    prev = 1;
#    prvek = srt[1];
#    vsz = length(srt);
#    prvku = length(meanDistribution);
#    i = arraysUsed + 1;
#    #print(i)
#    if(vsz>1)
#    {
#        for(j in 2:vsz)
#        {
#            pos = round((vsz-prvku+(prvku-1)*j)/(vsz-1))
#            if(pos==prev)
#            {
#                prvek=c(prvek, srt[j])
#            }else
#            {
#                meanDistribution[pos-1] = (1-1/i)*meanDistribution[pos-1] + mean(prvek)/i;
#                prvek = srt[j];
#                prev = pos;
#            }
#        }
#    }
#    meanDistribution[prev] = (1-1/i)*meanDistribution[prev] + mean(prvek)/i;
#    return(meanDistribution);
	prvku = length(meanDistribution)
	interpolatedVector = interpolateSortedVector(srt, prvku)
	i = arraysUsed + 1
	meanDistribution <- (1 - 1/i) * meanDistribution + (1/i) * interpolatedVector
	return(meanDistribution)
	
}

interpolateSortedVector <- function
###Interpolates given sorted vector to the vector of different
###length. It does not sort input vector thus for unsorted vectors do not guarantee functionality.
###Internal function.
##title<<Interpolate sorted vector
(vector,##<<Sorted vector to interpolate.
newSize##<<Size of the vector to produce.
)
{
	if (length(vector) == 0) {
		return(as.double(rep(NA, newSize)))
	}
	if (newSize == 0) {
		return(double())
	}
	#ret = .C("interpolateSortedVector_", as.double(vector), as.integer(length(vector)), 
	#		as.integer(newSize), newVector = double(newSize))$newVector
	ret = interpolateSortedVectorRcpp_(vector, newSize)
	return(ret)
}


meanDistribution <- function
###This function processes arrays in the object beadLevelData from package beadarray and returns sorted double vector.
###The vector has length prvku. And the distribution of this vector is a "mean" of all distributions of distributionChannel quantity in arrays.
###In case that probe numbers are different from prvku it does some averaging.
##title<<Produce sorted double vector with mean distribution.
(b,##<<Object beadLevelData from package beadarray or list of these objects
normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical 
##vector of the length equals to the number of arrays in b or list of such vectors 
##if b is a list of beadLevelData classes (defaults to NULL). 
distributionChannel="Grn",##<<Channel to do mean distribution from (defaults to "Grn").
channelInclude=NULL,##<<This field allows user to set channel with weights which have to be in {0,1}. 
##All zero weighted items are excluded from quantile normalization and the value asigned to such probes is a close to value which would be assigned to them if not being excluded. 
##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to NULL).
prvku##<<Number of items in a resulting double vector. Prvku must not be more than minimal number of indluded items in any distributionChannel.
)
{
    if(class(b)!="list")
    {
        b = list(b)
        if(!is.null(normalizationMod))
        {
            normalizationMod = list(normalizationMod)
        }
    }
    checkIntegrityOfListOfBeadLevelDataObjects(b, "warn")
    
    #channelInclude initialization
    if(is.null(channelInclude))
    {
        nowts=TRUE;
    }else
    {
        nowts=FALSE;
    }
    
    n <- numberOfDistributionElements(b, normalizationMod, channelInclude);
    if(missing(prvku) || is.null(prvku))
    {
        prvku <- n
    }
    if(prvku > n)
    {
        stop("Prvku must not be more than minimal number of included items in any distributionChannel.")
    }

    arraysUsed = 0;
    for(i in 1:length(b))
    {
        bbb = b[[i]];
        iteratorSet = 1:nrow(bbb);
        if(!is.null(normalizationMod))
        {
            iteratorSet = iteratorSet[normalizationMod[[i]]];
        }
        for(j in iteratorSet)
        {
            if(nowts)
            {
                srt = sort(bbb[[j]][,distributionChannel]);
            }else
            {
                srt = sort(bbb[[j]][bbb[[j]][,channelInclude]==1,][,distributionChannel])
            }
            
            if(arraysUsed>0)
            {
                meanDistribution <- updateMeanDistribution(meanDistribution, srt, arraysUsed)
            }else
            {
                meanDistribution <- initMeanDistribution(srt, prvku);
            }
            arraysUsed = arraysUsed + 1;
        }
    }

    if(is.unsorted(meanDistribution))
    {
        stop("UNSORTED SEQUENCE - PLEASE REPORT BUG!")
    }
    return(meanDistribution)
}

