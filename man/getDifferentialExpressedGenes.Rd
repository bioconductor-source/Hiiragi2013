\name{getDifferentialExpressedGenes}
\alias{getDifferentialExpressedGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find differentially expressed genes}
\description{
  Returns differentially expressed genes between two conditions.
}
\usage{
getDifferentialExpressedGenes(x, groups, g1, g2, theta = 0.5, FDRcutoff = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{ExpressionSet} containing a matrix of expression values with rows representing features and columns samples}
  \item{groups}{list of integer vectors specifying the grouping of samples}
  \item{g1}{character string specifying the name of the first tested group from \code{groups}}
  \item{g2}{character string specifying the name of the second tested group from \code{groups}}
  \item{theta}{numeric, probability with values in [0,1] used for quantile filtering of variance} 
  \item{FDRcutoff}{numeric, p-value cutoff}
}
\details{
  Differentailly expressed features are selected based on a t-test with the adjusted p-value cutoff specified by \code{FDRcutoff}.
  
  The filtering selects only features whose variance is greater than the quantile defined by the probability \code{theta}.
}
\value{Integer vector containing indices of differentially expressed features from \code{x}}
\author{Wolfgang Huber, 2013}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, groups, g1, g2, theta = 0.5, FDRcutoff = 0.05) 
{
    stopifnot(g1 \%in\% names(groups), g2 \%in\% names(groups))
    fac = rep(NA_integer_, ncol(x))
    fac[groups[[g1]]] = 1
    fac[groups[[g2]]] = 2
    fac = factor(fac)
    rv = rowVars(exprs(x)[, !is.na(fac)])
    passfilter = (rv > quantile(rv, probs = theta))
    tt = rowttests(x, fac)
    tt$adjp = rep(NA_real_, nrow(tt))
    tt$adjp[passfilter] = p.adjust(tt$p.value[passfilter])
    ord = order(tt$adjp)
    differentially = ord[tt$adjp[ord] <= FDRcutoff & !is.na(tt$adjp[ord])]
    return(differentially)
  }
}