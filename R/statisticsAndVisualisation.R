# Functions for bead level summary production. 
# Author: Vojtěch Kulvait
# Licence: GPL-3
###############################################################################

writeBackgroundImages <- structure(function
            ###This function writes images with background distribution according to foreground before and after background subtraction.
            ##title<<Write Background Images
            (b, ##<<Single beadLevelData object.
            spotsToGenerate=NULL,##<<NULL for generate images for all spots from b. Otherwise specifies logical vector of the length equals to the number of arrays in b with TRUE for images to generate. 
            imageType = c("jpg", "png", "eps"),##<<Type of images produced, either jpg, png or eps
            channelForeground="GrnF",##<<Name of channel of foreground.
            channelBackground="GrnB",##<<Name of channel of background.
            SDMultiple=3,##<<Correct on this level.
            includePearson=FALSE,##<<Include Pearson corelation.
            outputDir=getwd(),##<<Directory where to output images.
            width=505,##<<Width of image (default 505 fits well for 86mm 150dpi illustration in Bioinformatics journal:)
            height=505##<<Height of image
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn");
    singleChannelExistsIntegrityWithLogicalVector(b, spotsToGenerate, channelForeground, "error")
    singleChannelExistsIntegrityWithLogicalVector(b, spotsToGenerate, channelBackground, "error")
    iteratorSet = 1:nrow(b)
    imageType = match.arg(imageType)
    if(!is.null(spotsToGenerate))
    {
        singleCheckIntegrityLogicalVector(spotsToGenerate, b, "error")
        iteratorSet = iteratorSet[spotsToGenerate];
    }
    for(i in iteratorSet)
    {
        if(imageType == "jpg")
        {
            jpeg(paste(file.path(outputDir, as.character(b@sectionData$SampleGroup[i,])), ".jpg", sep=""), width = width, height = height, quality = 100)
        }else if(imageType == "png")
        {
            png(paste(file.path(outputDir, as.character(b@sectionData$SampleGroup[i,])), ".png", sep=""), width = width, height = height)
        }else if(imageType == "eps")
        {
            setEPS()
            postscript(paste(file.path(outputDir, as.character(b@sectionData$SampleGroup[i,])), ".eps", sep=""), horizontal = FALSE, onefile = FALSE)
        }
        plotBackgroundImageBeforeCorrection(b, i, channelForeground, channelBackground, includePearson)
        dev.off()
        
        if(imageType == "jpg")
        {
            jpeg(paste(file.path(outputDir,paste(as.character(b@sectionData$SampleGroup[i,]),"_CORRECTED", sep="")), ".jpg", sep=""), width = width, height = height, quality = 100)
        }else if(imageType == "png")
        {
            png(paste(file.path(outputDir,paste(as.character(b@sectionData$SampleGroup[i,]),"_CORRECTED", sep="")), ".png", sep=""), width = width, height = height)
        }else if(imageType == "eps")
        {
            setEPS()
            postscript(paste(file.path(outputDir,paste(as.character(b@sectionData$SampleGroup[i,]),"_CORRECTED", sep="")), ".eps", sep=""), horizontal = FALSE, onefile = FALSE)
        }
        plotBackgroundImageAfterCorrection(b, i, channelForeground, channelBackground, SDMultiple, includePearson)
        dev.off()
    }
}, ex=function()
{
    if(require("blimaTestingData") && interactive())
    {
        #Write background images before and after correction for background into /tmp directory. This function creates two jpg images for condition D. Output files are /tmp/6898481102_D_CORRECTED.jpg and /tmp/6898481102_D.jpg.
        data(blimatesting)
        p = pData(blimatesting[[2]]@experimentData$phenoData)
        spotsToGenerate = p$Group %in% "D";
        writeBackgroundImages(blimatesting[[2]], imageType="jpg", spotsToGenerate=spotsToGenerate, includePearson=FALSE, outputDir="/tmp", width=505, height=505)
    }else
    {
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})

plotBackgroundImageBeforeCorrection <- structure(function
        ###This function plots image of background distribution versus to foreground before background subtraction.
        ##title<<Plot background image before correction
        (b, ##<<Single beadLevelData object.
        index,##<<Index of spot to generate.
        channelForeground="GrnF",##<<Name of channel of foreground.
        channelBackground="GrnB",##<<Name of channel of background.
        includePearson=FALSE##<<Include Pearson corelation.
    )
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn");
    if(missing(index) || !is.numeric(index))
    {
        stop("Integer field index is mandatory.")
    }
    if(includePearson)
    {
        pearson <- cor(b[[index]][,channelBackground], b[[index]][, channelForeground])
        title <- paste("Background vs foreground, pearson cor=", round(pearson, 3), sep="")
    }else
    {
        title <- "Background vs foreground"
    }
    graphics::plot(b[[index]][,channelBackground], b[[index]][, channelForeground], xlab="Background", ylab="Foreground", main=title, pch=21)
}, ex=function()
{
    if(require("blimaTestingData") && interactive())
    {
        #Write background images before correction. This function prints graph for condition D4. Call dev.off() to close.
        data(blimatesting)
        p = pData(blimatesting[[2]]@experimentData$phenoData)
        index = base::match("D4", p$Name)
        plotBackgroundImageBeforeCorrection(blimatesting[[2]], index)
    }else
    {
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})

plotBackgroundImageAfterCorrection <- structure(function
        ###This function plots image of background distribution versus to foreground after background subtraction.
        ##title<<Plot background image after correction
        (b, ##<<Single beadLevelData object.
        index,##<<Index of spot to generate.
        channelForeground="GrnF",##<<Name of channel of foreground.
        channelBackground="GrnB",##<<Name of channel of background.
        SDMultiple=3,##<<Correct on this level.
        includePearson=FALSE##<<Include Pearson corelation.
    )
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn");
    if(missing(index) || !is.numeric(index))
    {
        stop("Integer field index is mandatory.")
    }
    men <- mean(b[[index]][,channelBackground])
    sd <- sd(b[[index]][,channelBackground])
    correctedVEC <-filterBg(b[[index]][,channelBackground], SDMultiple)
    lg1 <- sum(correctedVEC)
    lg2 <- length(correctedVEC)
    pomer <- (lg1/lg2)
    
    if(includePearson)
    {
        pearson <- cor(b[[index]][,channelBackground][correctedVEC], b[[index]][, channelForeground][correctedVEC])
        title <- paste("Background corrected, remains ",round(pomer*100, 2),"%, pearson cor=", round(pearson, 3), sep="")
    }else
    {
        title <- paste("Background corrected, remains ",round(pomer*100, 2),"%.", sep="")
    }
    graphics::plot(b[[index]][,channelBackground][correctedVEC], b[[index]][, channelForeground][correctedVEC], xlab="Background", ylab="Foreground", main=title, pch=21)
}, ex=function()
{
    if(require("blimaTestingData") && interactive())
    {
        #Write background images after correction. This function prints graph for condition D4. Call dev.off() to close.
        data(blimatesting)
        p = pData(blimatesting[[2]]@experimentData$phenoData)
        index = base::match("D4", p$Name)
        plotBackgroundImageAfterCorrection(blimatesting[[2]], index)
    }else
    {
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})
        
chipArrayStatistics <- structure(function
        ###This function returns table with statistics of single beadLevelData object indexed by order of spots.
        ###It prints number of beads on each array spot mean foreground intensity and optionally mean background intensity, mean number of 
        ###beads in probe set and unbiased estimate of standard deviations of these parameters.
        ###Optionaly you can also obtain percentage of removed beads within excludedOnSDMultiple
        ###multiple of standard deviations from the background value.
        ##title<<Statistics of beadLevelData
        (b,##<<Single beadLevelData object.
        includeBeadStatistic = TRUE,##<<Include number of beads per probe in output.
        channelForeground="GrnF",##<<Name of channel of foreground.
        channelBackground="GrnB",##<<Name of channel of background.
        includeBackground=TRUE,##<<Whether to output background data.
        excludedOnSDMultiple=NA##<<If positive number, print how much percents of the background lies more than excludedOnSDMultiple multipliers of standard deviation estimate away from background mean.
)
{
    checkIntegrityOfSingleBeadLevelDataObject(b, "warn");
    singleChannelExistsIntegrityWithLogicalVector(b, NULL, channelForeground, "error")
    singleChannelExistsIntegrityWithLogicalVector(b, NULL, channelBackground, "error")
    dsc <- matrix(0, nrow(b), 0);
    beads <- numeric();
    for(i in 1:nrow(b))
    {
        beads = c(beads, nrow(b[[i]]));
    }
    dsc <- insertColumn(dsc, beads, "Beads")
    
    fg <- numeric();
    fgs <- numeric();
    for(i in 1:nrow(b))
    {
        fg = c(fg, mean(b[[i]][,channelForeground]));
        fgs = c(fgs,  sd(b[[i]][,channelForeground]))
    }    
    dsc <- insertColumn(dsc, fg, "Mean FG")
    dsc <- insertColumn(dsc, fgs, "SD FG")
    
    if(includeBackground)
    {
        bg <- numeric();
        bgs <- numeric();
        for(i in 1:nrow(b))
        {
            bg = c(bg, mean(b[[i]][,channelBackground]));
            bgs = c(bgs,  sd(b[[i]][,channelBackground]))
        }    
        dsc <- insertColumn(dsc, bg, "Mean BG")
        dsc <- insertColumn(dsc, bgs, "SD BG")
        
        if(!is.na(excludedOnSDMultiple) && excludedOnSDMultiple>0)
        {
            pct <- numeric();
            for(i in 1:nrow(b))
            {
                cvn = filterBg(b[[i]][,channelBackground], excludedOnSDMultiple)
                pct = c(pct, (length(cvn)-sum(cvn))*100/length(cvn));
            }    
        }
    }
    
    if(includeBeadStatistic)
    {
        bg <- numeric();
        bgs <- numeric();
        for(i in 1:nrow(b))
        {
            bg = c(bg, mean(table(b[[i]][,"ProbeID"])));
            bgs = c(bgs, sd(table(b[[i]][,"ProbeID"])));
        }    
        dsc <- insertColumn(dsc, bg, "Mean BPP")
        dsc <- insertColumn(dsc, bgs, "SD BPP")
    }
    
    if(!is.na(excludedOnSDMultiple) && excludedOnSDMultiple>0)
    {
        dsc <- insertColumn(dsc, pct, "PCT")
    }
    return(dsc);
}, ex=function()
{
    if(require("blimaTestingData") && interactive())
    {
        #To print basic statistic data about blimatesting[[1]] object.
        data(blimatesting)
        array1stats = chipArrayStatistics(blimatesting[[1]], includeBeadStatistic=TRUE,
                excludedOnSDMultiple=3)
        array1pheno = pData(blimatesting[[1]]@experimentData$phenoData)
        array1stats = data.frame(array1pheno$Name, array1stats)
        colnames(array1stats)[1] <- "Array";
        print(array1stats);
    }else
    {
        print("To run this example, please install blimaTestingData package from bioconductor by running BiocManager::install('blimaTestingData').");
    }
})


insertColumn <- function
###Internal
##title<<Internal function to support chipArrayStatistics
(matrix,##<<Object to insert column to
column,##<<Column to insert
name##<<Name of column to assign.
)
{
    matrix <- cbind(matrix, column);
    colnames(matrix)[ncol(matrix)] <- name;
    return(matrix);
}
