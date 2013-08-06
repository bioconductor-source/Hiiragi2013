\name{myHeatmap}
\alias{myHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot heatmaps
}
\description{
Plots a heatmap.
}
\usage{
myHeatmap(x, collapseDuplicateFeatures = FALSE, haveColDend = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ExpressionSet}
  \item{collapseDuplicateFeatures}{logical, if \code{TRUE} multiple features per gene will be averaged}
  \item{haveColDend}{logical}
}
\value{
The result of the function is a plot.
}
\author{Wolfgang Huber, 2013}
\seealso{\code{\link{myHeatmap2}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--  or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, collapseDuplicateFeatures = FALSE, haveColDend = FALSE) 
{
    mat = scale(t(exprs(x)))
    colnames(mat) = fData(x)[, "symbol"]
    rownames(mat) = sub(".CEL$", "", rownames(mat))
    if (collapseDuplicateFeatures) {
        sp = split(seq(along = colnames(mat)), f = colnames(mat))
        mat = sapply(sp, function(j) rowMeans(mat[, j, drop = FALSE]))
    }
    rowdd = as.dendrogram(hclust(dist(mat)))
    roword = order.dendrogram(rowdd)
    coldd = as.dendrogram(hclust(dist(t(mat))))
    colord = order.dendrogram(coldd)
    print(levelplot(mat[roword, colord], aspect = "fill", xlab = "", 
        ylab = "", scales = list(x = list(rot = 90), raster = TRUE), 
        col.regions = bluered(100), colorkey = list(space = "left", 
            height = 0.15, useRaster = TRUE), legend = if (haveColDend) 
            list(top = list(fun = dendrogramGrob, args = list(x = rowdd, 
                side = "top", size = 10)))))
  }
}