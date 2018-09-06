# Functions for bead level summary production. 
# Author: Vojtěch Kulvait
# Licence: GPL-3
###############################################################################

varianceBeadStabilise <- structure(function
        ###This function does variance stabilising step on bead level.
        ##title<<Bead level VST.
        (b,##<<List of beadLevelData objects (or single object).
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equal to the number of arrays in b or list of such vectors if b is a list of beadLevelData classes. 
        quality="qua",##<<Quality to analyze, default is "qua".
        channelInclude="bgf",##<<This field allows user to set channel with weights which have to be in {0,1}.
        ##All zero weighted items are excluded from t-test.
        ##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to "bgf").
        channelOutput="vst"##<<Output from VST.
)
{
    #For vst first you have to provide vector of sumarized values together with their standard deviations
    waslist = checkIntegrity(b, "warn")
    if(!waslist)
    {
        b = list(b)
        if(!is.null(normalizationMod))
        {
            normalizationMod = list(normalizationMod)
        }
    }
    channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, quality, "error")
    if(!is.null(channelInclude))
    {
        channelExistsIntegrityWithLogicalVectorList(b, normalizationMod, channelInclude, "error")
    }
    x=list();
    for(i in 1:length(b))
    {
        x=c(x, varianceBeadStabiliseSingleArray(b[[i]], normalizationMod[[i]], quality, channelInclude, channelOutput));
    }
    return(x);
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
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})

varianceBeadStabiliseSingleArray <- function
        ###This function is not intended to direct use it takes single beadLevelData object and do bead level variance stabilisation.
        ##title<<Bead level VST.
        (b,##<<Object beadLevelData.
        normalizationMod=NULL,##<<NULL for normalization of all input b. Otherwise specifies logical vector of the length equals to the number of arrays in b. 
        quality="qua",##<<Quality to analyze, default is "qua".
        channelInclude="bgf",##<<This field allows user to set channel with weights which have to be in {0,1}. 
        ##All zero weighted items are excluded from t-test. 
        ##You can turn this off by setting this NULL. This option may be used together with bacgroundCorrect method or/and with beadarray QC (defaults to "bgf").
        channelOutput="vst"##<<Output from VST.
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn")
    if(is.null(channelInclude))
    {
        nowts=TRUE;
    }else
    {
        nowts=FALSE;
    }
    iteratorSet = 1:nrow(b);
    if(!is.null(normalizationMod))
    {
        iteratorSet = iteratorSet[normalizationMod];
    }
    for(i in iteratorSet)
    {
        if(nowts)
        {
            meansdtable <- aggregateAndPreprocess(b[[i]][,c("ProbeID", quality)], quality)
        }else
        {
            meansdtable <- aggregateAndPreprocess(b[[i]][b[[i]][,channelInclude]==1,c("ProbeID", quality)], quality)
        }
        c3estimate = median(meansdtable[order(meansdtable[,"mean"]),][1:nrow(meansdtable),"sd"], na.rm=TRUE)^2
        #this is hack due to incorrect implementation of vst to pass variance as sd
        vstresults <- vstFromLumi(meansdtable[,"mean"], meansdtable[,"sd"], backgroundStd=c3estimate)
        parameters <- attr(vstresults, "parameter")
        fun <- attr(vstresults, "transformFun")
        src = b[[i]][,quality]
        if(fun=="log")
        {
            todst <- parameters[["g"]] * log(parameters[["a"]] + parameters[["b"]]*src) + parameters[["Intercept"]]
        }else
        {
            todst <- parameters[["g"]] * asinh(parameters[["a"]] + parameters[["b"]]*src) + parameters[["Intercept"]]
        }
        b = setWeights(b, wts=todst, array=i, wtName=channelOutput)
    }
    return(b)
}


vstFromLumi <- function
###This function is derived from copy and paste of lumi::vst function. Since lumi package has extensive imports I decided to hardcode this function to the blima instead of importing lumi package.
##title<<Function from LGPL lumi package 2.16.0 
(u,##<<The mean of probe beads
std,##<<The standard deviation of the probe beads
nSupport=min(length(u), 500),##<<Something for c3 guess.
backgroundStd=NULL,##<<Estimate the background variance c3. Input should be variance according to article, not SD.
lowCutoff=1/3##<<Something for c3 guess.
)
##references<< \url{http://www.bioconductor.org/packages/release/bioc/html/lumi.html}
##author<<authors are Pan Du, Simon Lin, the function was edited by
{
	# u is the mean of probe beads
	# std is the standard deviation of the probe beads
	## Estimate the background variance c3
	c3 <- ifelse (is.null(backgroundStd), 0, backgroundStd)
	
	ord <- order(u); u.bak <- u
	u <- u[ord]; std <- std[ord]
	
	## remove NAs if exists
	na.ind <- which(is.na(u) | is.na(std))
	if (length(na.ind) > 0) {
		u <- u[-na.ind]
		std <- std[-na.ind]
	}
	
	if (any(std < 0)) {
		stop('Negative expression standard deviation is not allowed!')
	}
	
	
	
    ## downsampling to speed up 
    # if (min(u) < 1) {
    #   offset <- 1 - min(u)
    # } else {
    #   offset <- 0
    # }
    offset <- 1 - min(u)
    downSampledU <- 2^seq(from=log2(min(u + offset)), to=log2(max(u + offset)), length=nSupport) - offset
    
    minU <- log2(max(100 + offset, min(u)))
    maxU <- log2(max(u))
    # uCutoff <- 2^((maxU + minU)/2)
    uCutoffLow <- 2^(minU + (maxU - minU) * lowCutoff)
    uCutoffHigh <- 2^(minU + (maxU - minU) * 4/5)
    selInd <- (u > uCutoffLow & u < uCutoffHigh)
    selLowInd <- (u < uCutoffLow)
    if (c3 != 0) {
    	selInd <- selInd & (std^2 > c3)
    	dd <- data.frame(y=sqrt(std[selInd]^2 - c3), x1=u[selInd])
    	# if (nrow(dd) > 5000) dd <- dd[sample(1:nrow(dd), 5000),]
    	lmm <- lm(y ~ x1, dd)
    	c1 <- lmm$coef[2]
    	c2 <- lmm$coef[1]
    } else {
    	## use iteraction to estimate the parameters when background level is not known
    	iterNum <- 0
    	c3.i <- 0
    	while(iterNum <= 20) {
    		selInd.i <- selInd & (std^2 > c3.i)
    		dd <- data.frame(y=sqrt(std[selInd.i]^2 - c3.i), x1=u[selInd.i])
    		# if (nrow(dd) > 5000) dd <- dd[sample(1:nrow(dd), 5000),]
    		lm.i <- lm(y ~ x1, dd)
    		c1.i <- lm.i$coef[2]
    		c2.i <- lm.i$coef[1]
    		y <- std[selLowInd]
    		x <- u[selLowInd]
    		cc <- y^2 - (c1.i * x + c2.i)^2
    		c3.i.new <- mean(cc, trim=0.05)
    		if (c3.i.new < 0) {
    			break
    		} else {
    			if (abs(c3.i.new - c3.i) < 1e-5) break
    			c3.i <- c3.i.new
    		}
    		iterNum <- iterNum + 1
    	}
    	c1 <- c1.i; c2 <- c2.i; c3 <- c3.i
    	if (c3 < 0) c3 <- 0
    }
    smoothStd <- ((c1 * downSampledU + c2)^2 + c3)^(1/2)
	
	## calculate the integration (h function is the integral)
	
    if (c3 == 0) {
    	## Transform function h(x) = g * log(a + b * x)
    	g <- 1/c1
    	a <- c2
    	b <- c1
    	tmp <- a + b * u.bak
    	if (any(tmp < 0)) {
    		transformedU <- log(u.bak)
    		g <- 1; a <- 0; b <- 1
    	} else {
    		transformedU <- g * log(a + b * u.bak)
    	}
    	transFun <- 'log'
    } else {
    	## Transform function h(x) = g * asinh(a + b * x)
    	g <- 1/c1
    	a <- c2/sqrt(c3)
    	b <- c1/sqrt(c3)
    	transformedU <- g * asinh(a + b * u.bak)
    	transFun <- 'asinh'
    }
    transform.parameter <- c(a, b, g, 0)
    names(transform.parameter) <- c('a', 'b', 'g', 'Intercept')
	
	
	cutInd <- which.min(abs(u.bak - uCutoffLow))
	maxInd <- which.max(u.bak)
	y <- c(u.bak[cutInd], u.bak[maxInd])
	x <- c(transformedU[cutInd], transformedU[maxInd])
	m <- lm(log2(y) ~ x)
	
    transform.parameter <- c(a, b, g * m$coef[2], m$coef[1])
    names(transform.parameter) <- c('a', 'b', 'g', 'Intercept')
    ## The transform parameter is in the transFun below
    ## transFun <- g * asinh(a + b * x) * m$coef[2] + m$coef[1]
	
	transformedU <- predict(m, data.frame(x=transformedU))
	attr(transformedU, 'parameter') <- transform.parameter
	attr(transformedU, 'transformFun') <- transFun
	
	return(transformedU)
}


