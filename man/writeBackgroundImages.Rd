\name{writeBackgroundImages}
\alias{writeBackgroundImages}
\title{Write Background Images}
\description{This function writes images with background distribution according to foreground before and after background subtraction.}
\usage{writeBackgroundImages(b, spotsToGenerate = NULL, imageType = c("jpg", 
    "png", "eps"), channelForeground = "GrnF", channelBackground = "GrnB", 
    SDMultiple = 3, includePearson = FALSE, outputDir = getwd(), 
    width = 505, height = 505)}
\arguments{
  \item{b}{Single beadLevelData object.}
  \item{spotsToGenerate}{NULL for generate images for all spots from b. Otherwise specifies logical vector of the length equals to the number of arrays in b with TRUE for images to generate.}
  \item{imageType}{Type of images produced, either jpg, png or eps}
  \item{channelForeground}{Name of channel of foreground.}
  \item{channelBackground}{Name of channel of background.}
  \item{SDMultiple}{Correct on this level.}
  \item{includePearson}{Include Pearson corelation.}
  \item{outputDir}{Directory where to output images.}
  \item{width}{Width of image (default 505 fits well for 86mm 150dpi illustration in Bioinformatics journal:)}
  \item{height}{Height of image}
}



\author{Vojtech Kulvait}




\examples{
if(require("blimaTestingData") && interactive())
{
    #Write background images before and after correction for background into /tmp directory. This function creates two jpg images for condition D. Output files are /tmp/6898481102_D_CORRECTED.jpg and /tmp/6898481102_D.jpg.
    data(blimatesting)
    p = pData(blimatesting[[2]]@experimentData$phenoData)
    spotsToGenerate = p$Group \%in\% "D";
    writeBackgroundImages(blimatesting[[2]], imageType="jpg", spotsToGenerate=spotsToGenerate, includePearson=FALSE, outputDir="/tmp", width=505, height=505)
}else
{
    print("To run this example, please install blimaTestingData package from bioconductor by running biocLite('blimaTestingData').");
}
}
