\name{myHeatmap}
\alias{myHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Heatmaps}
\description{
Performs hierarchical cluster analysis of features and samples from an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} and plots the resulting heatmap.
}
\usage{
myHeatmap(x, collapseDuplicateFeatures = FALSE, haveColDend = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing a matrix of expression values with rows representing features and columns samples}
  \item{collapseDuplicateFeatures}{logical, if \code{TRUE} multiple features per gene will be averaged}
  \item{haveColDend}{logical, if \code{TRUE} displays the grouping of samples using a dendrogram}
}
\value{
The result of the function is a plot.
}
\author{Wolfgang Huber, 2013}
\seealso{\code{\link{myHeatmap2}}}
\examples{
data(x)
groups = with(pData(x), list(
  "E3.25"       = which(genotype=="WT" & Embryonic.day=="E3.25"),
  "E3.5 (EPI)"  = which(genotype=="WT" & Embryonic.day=="E3.5" & lineage=="EPI"),
  "E3.5 (PE)"   = which(genotype=="WT" & Embryonic.day=="E3.5" & lineage=="PE")))
samples = unlist(groups)

# heatmap of differentially expressed genes along the transition from E3.25 to E3.5
de = union(
  getDifferentialExpressedGenes(x, groups, "E3.25", "E3.5 (EPI)"),
  getDifferentialExpressedGenes(x, groups, "E3.25", "E3.5 (PE)"))
myHeatmap(x[de, samples], collapseDuplicateFeatures=TRUE)
}