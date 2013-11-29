\name{getDifferentialExpressedGenes}
\alias{getDifferentialExpressedGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Differentially Expressed Genes}
\description{
  Returns differentially expressed genes between two conditions.
}
\usage{
getDifferentialExpressedGenes(x, groups, g1, g2, theta = 0.5, FDRcutoff = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing a matrix of expression values with rows representing features and columns samples}
  \item{groups}{list of integer vectors specifying the grouping of samples}
  \item{g1}{character string specifying the name of the first tested group from \code{groups}}
  \item{g2}{character string specifying the name of the second tested group from \code{groups}}
  \item{theta}{numeric, probability with values in [0,1] used for quantile filtering of variance} 
  \item{FDRcutoff}{numeric, p-value cutoff}
}
\details{
  Differentially expressed features are selected based on a t-test with the adjusted p-value cutoff specified by \code{FDRcutoff}.
  
  The filtering selects only features whose variance is greater than the quantile defined by the probability \code{theta}.
}
\value{Integer vector containing indices of differentially expressed features from \code{x}.}
\author{Wolfgang Huber, 2013}
\examples{
data(x)
groups = with(pData(x), list(
  "E3.25"       = which(genotype=="WT" & Embryonic.day=="E3.25"),
  "E3.5 (EPI)"  = which(genotype=="WT" & Embryonic.day=="E3.5" & lineage=="EPI"),
  "E3.5 (PE)"   = which(genotype=="WT" & Embryonic.day=="E3.5" & lineage=="PE")))
  
# get a list of differentially expressed genes along the transition from E3.25 to E3.5
de = union(
  getDifferentialExpressedGenes(x, groups, "E3.25", "E3.5 (EPI)"),
  getDifferentialExpressedGenes(x, groups, "E3.25", "E3.5 (PE)"))
fData(x[de,])$symbol
}