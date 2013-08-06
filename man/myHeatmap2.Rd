\name{myHeatmap2}
\alias{myHeatmap2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot heatmaps
}
\description{
Plots a heatmap with clustering of rows and colums specified by the \code{rowGroups} and \code{colGroups}, respectively.
}
\usage{
myHeatmap2(x, rowGroups = factor(rep(1, nrow(x))), colGroups = factor(rep(1, ncol(x))), keeprownames = TRUE, colors = colorRampPalette(brewer.pal(9, "Blues")[-1])(100), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ExpressionSet}
  \item{rowGroups}{factor which length aligns with rows of x}
  \item{colGroups}{factor which length aligns with columns of x}
  \item{keeprownames}{logical, if true displays row names}
  \item{colors}{color palette used to plot the heatmap}
  \item{\dots}{arguments passed to the internal \code{ordermat} function}
}
\value{
The result of the function is a plot.
}
\author{Wolfgang Huber, 2013}
\seealso{\code{\link{myHeatmap}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, rowGroups = factor(rep(1, nrow(x))), colGroups = factor(rep(1, 
    ncol(x))), keeprownames = TRUE, colors = colorRampPalette(brewer.pal(9, 
    "Blues")[-1])(100), ...) 
{
    x = ordermat(x, rowGroups, ...)
    x = ordermat(t(x), colGroups, ...)
    if (!keeprownames) 
        colnames(x) = NULL
    print(levelplot(x, aspect = "fill", xlab = "", ylab = "", 
        scales = list(x = list(rot = 90), raster = TRUE), col.regions = colors, 
        colorkey = list(space = "left", height = 0.15, useRaster = TRUE)))
  }
}
