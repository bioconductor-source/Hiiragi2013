pamCluster <-
function(ngenes, x, k=2) {
  topGenes = order(rowVars(exprs(x)), decreasing=TRUE)[seq_len(ngenes)]
  dists = dist(t(exprs(x)[topGenes, ]))
  pam(dists, k=k, stand=FALSE, cluster.only=TRUE)
}
