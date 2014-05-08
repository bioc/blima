\name{chipArrayStatistics}
\alias{chipArrayStatistics}
\title{Statistics of beadLevelData}
\description{This function returns table with statistics of single beadLevelData object indexed by order of spots.
It prints number of beads on each array spot mean foreground intensity and optionally mean background intensity, mean number of 
beads in probe set and unbiased estimate of standard deviations of these parameters.
Optionaly you can also obtain percentage of removed beads within excludedOnSDMultiple
multiple of standard deviations from the background value.}
\usage{chipArrayStatistics(b, includeBeadStatistic = TRUE, channelForeground = "GrnF", 
    channelBackground = "GrnB", includeBackground = TRUE, excludedOnSDMultiple = NA)}
\arguments{
  \item{b}{Single beadLevelData object.}
  \item{includeBeadStatistic}{Include number of beads per probe in output.}
  \item{channelForeground}{Name of channel of foreground.}
  \item{channelBackground}{Name of channel of background.}
  \item{includeBackground}{Whether to output background data.}
  \item{excludedOnSDMultiple}{If positive number, print how much percents of the background lies more than excludedOnSDMultiple multipliers of standard deviation estimate away from background mean.}
}



\author{Vojtech Kulvait}




\examples{
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
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
}
}
