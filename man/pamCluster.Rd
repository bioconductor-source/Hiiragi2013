\name{pamCluster}
\alias{pamCluster}
\title{Clustering of Most Variable Genes}
\description{
The function \code{pamCluster} selects the \code{ngenes} most variable genes and performs their clustering using the partitioning around medoids method \code{\link[cluster:pam]{pam}}.
}
\usage{
pamCluster(ngenes, x, k = 2)
}
\arguments{
  \item{ngenes}{numeric, the number of most variable genes to select} 
  \item{x}{\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} containing gene expression values}
  \item{k}{positive integer specifying the number of clusters}
}
\value{
  Integer vector specifying the clustering.
}
\seealso{pam}
\author{Wolfgang Huber}
\examples{
data("x")
y = x[, x$Embryonic.day=="E3.5"]

## perform the clustering 
pc = pamCluster(50, y, k=3)

## display clustering vs. sample lineage
plot(as.factor(pData(y)$lineage), pc, yaxt="n", xlab="lineage", ylab="cluster")
}
