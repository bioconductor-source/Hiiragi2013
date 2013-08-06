\name{pamCluster}
\alias{pamCluster}
\title{
Clustering of the most variable genes
}
\description{
The function \code{pamCluster} selects the \code{ngenes} most variable genes in
the data matrix \code{x} and seeks for two clusters using the
partitioning around medoids method \code{pam}.
}
\usage{
pamCluster(ngenes, x, k = 2)
}
\arguments{
  \item{ngenes}{numeric, the number of most variable genes to select} 
  \item{x}{\code{ExpressionSet} containing gene expression values}
  \item{k}{positive integer specifying the number of clusters}
}
\value{
  Integer vector specifying the clustering.
}
\seealso{pam}
\author{Wolfgang Huber}
\examples{
  data("x")
  pc = pamCluster(40, x, k=4)
}
