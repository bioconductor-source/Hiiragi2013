\name{plotProjection}
\alias{plotProjection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting overview of the sample expression}
\description{
Projection of sample expression profiles on
the differential expression signature.
}
\usage{
plotProjection(projection, label, col, colourMap)
}
\arguments{
  \item{projection}{
%%     ~~Describe \code{projection} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{colourMap}{
%%     ~~Describe \code{colourMap} here~~
}
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
function (projection, label, col, colourMap) 
{
    projRank = (rank(projection) - 1)/(length(projection) - 1)
    gp = gpar(col = col)
    yticks = pretty(projection)
    yscale = yticks[c(1, length(yticks))]
    pushViewport(viewport())
    pushViewport(viewport(yscale = yscale, xscale = c(-0.25, 
        1), x = 0.4, width = 0.75, height = 0.95, clip = "off"))
    grid.lines(x = unit(c(0, 0), "native"), y = unit(yscale, 
        "native"))
    grid.points(x = unit(rep(0, length(projection)), "native"), 
        y = unit(projection, "native"), pch = 16, size = unit(0.75, 
            "char"), gp = gp)
    grid.text(x = unit(0.33, "native"), y = unit(projRank, "npc"), 
        label = label, gp = gp, just = c("left", "center"))
    grid.text(x = unit(-0.1, "native"), y = unit(yticks, "native"), 
        label = paste(yticks), hjust = 1)
    grid.segments(x0 = unit(0, "native"), x1 = unit(0.33, "native"), 
        y0 = unit(projection, "native"), y1 = unit(projRank, 
            "npc"), gp = gp)
    popViewport()
    pushViewport(viewport(yscale = c(0, length(colourMap) + 1), 
        xscale = c(0, 1), x = 0.85, width = 0.29, height = 0.2))
    grid.text(x = unit(0.3, "native"), y = unit(seq(along = colourMap), 
        "native"), names(colourMap), just = "left")
    grid.rect(x = unit(0.13, "native"), y = unit(seq(along = colourMap), 
        "native"), width = unit(0.26, "native"), height = unit(0.8, 
        "native"), gp = gpar(fill = colourMap))
    popViewport(2)
  }
}
