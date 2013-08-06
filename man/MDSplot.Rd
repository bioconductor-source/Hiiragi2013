\name{MDSplot}
\alias{MDSplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Scaling Plots
}
\description{
Visualization of data similarity using non-metric multidimensional scaling.
}
\usage{
MDSplot(x, mask, flip = integer(0), rotation = 0, cex = 2, col = x$sampleColour, panellabel, pointlabel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{ExpressionSet} containing a matrix of expression values with rows representing features and columns samples}
  \item{mask}{logical vector specifying the masking of datapoints; if missing all points are drawn}
  \item{flip}{integer vector containing indices of columns of the distance matrix whose sign should be flipped}
  \item{rotation}{numeric, angle in radians by which the plot should be rotated}
  \item{cex}{numeric, the value giving the amount by which text and symbols should be scaled relative to the current setting}
  \item{col}{character vector specifing the colors of points}
  \item{panellabel}{character, the label of the figure panel}
  \item{pointlabel}{character vector, the labels corresponding to samples; if present labels are drawn instead of points}
}
\value{
The result of the function is a plot.
}
\author{Wolfgang Huber, 2013}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, mask, flip = integer(0), rotation = 0, cex = 2, 
    col = x$sampleColour, panellabel, pointlabel) 
{
    mat = exprs(x)
    distances = dist(t(mat))
    mds = isoMDS(distances)$points
    for (i in flip) mds[, i] = -mds[, i]
    mds = mds \%*\% rotmat(rotation)
    par(mai = c(0.1, 0.1, 0.1, 0.1))
    xlim = range(mds[, 1])
    ylim = range(mds[, 2])
    if (missing(mask)) 
        mask = rep(TRUE, ncol(x))
    theColours = paste(col[mask], "a0", sep = "")
    plot(mds[mask, ], pch = 16, cex = cex, xlab = "", ylab = "", 
        xaxt = "n", yaxt = "n", asp = 1, xlim = xlim, ylim = ylim, 
        col = if (missing(pointlabel)) 
            theColours
        else par("bg"))
    if (missing(pointlabel)) {
        points(mds[mask, ], pch = 1, cex = cex, col = grey(0.3))
    }
    else {
        text(mds[mask, ], labels = pointlabel, adj = c(0.5, 0.5), 
            cex = cex * 2/3, col = theColours)
    }
    if (!missing(panellabel)) 
        text(par("usr")[1], par("usr")[4], panellabel, adj = c(-1, 
            1.5), cex = cex * 2/3, col = "#808080")
  }
}