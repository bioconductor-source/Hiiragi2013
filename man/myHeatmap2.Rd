\name{myHeatmap2}
\alias{myHeatmap2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Heatmaps}
\description{
Plots a heatmap with clustering of rows and columns specified by the \code{rowGroups} and \code{colGroups}, respectively.
}
\usage{
myHeatmap2(x, rowGroups = factor(rep(1, nrow(x))), colGroups = factor(rep(1, ncol(x))), keeprownames = TRUE, colors = colorRampPalette(brewer.pal(9, "Blues")[-1])(100), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing a matrix of expression values with rows representing features and columns samples}
  \item{rowGroups}{factor which length aligns with rows of \code{x}}
  \item{colGroups}{factor which length aligns with columns of \code{x}}
  \item{keeprownames}{logical, if \code{TRUE} displays row names}
  \item{colors}{color palette used to plot the heatmap}
  \item{\dots}{arguments passed to the internal \code{ordermat} function}
}
\value{
The result of the function is a plot.
}
\author{Wolfgang Huber, 2013}
\seealso{\code{\link{myHeatmap}}}
\examples{
## The function is currently defined as
function (x, rowGroups = factor(rep(1, nrow(x))), colGroups = factor(rep(1, ncol(x))),
    keeprownames = TRUE, colors = colorRampPalette(brewer.pal(9, "Blues")[-1])(100), ...) {
    x = ordermat(x, rowGroups, ...)
    x = ordermat(t(x), colGroups, ...)
    if (!keeprownames) 
        colnames(x) = NULL
    print(levelplot(x, aspect = "fill", xlab = "", ylab = "", 
        scales = list(x = list(rot = 90), raster = TRUE), col.regions = colors, 
        colorkey = list(space = "left", height = 0.15, useRaster = TRUE)))
  }
}
