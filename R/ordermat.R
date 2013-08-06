ordermat <-
function(x, fac, ordering="dend") {
  stopifnot(length(fac)==nrow(x), is.factor(fac))
  sp = split(seq(along=fac), fac)
  nx = lapply(seq(along=sp), function(i) {
    xi = x[sp[[i]],,drop=FALSE]
    if(nrow(xi)>1) {
      xi = switch(ordering,
        "none" = xi,
        "dend" = xi[order.dendrogram(as.dendrogram(hclust(dist(xi)))),,drop=FALSE],
        stop(sprintf("Unknown ordering method '%s'.", ordering)))
    }
    return(xi)
  })
  ngap = ceiling(sum(sapply(nx, nrow))/200)
  for(i in seq(along=nx)[-1])
    nx[[i]] = rbind(matrix(NA_real_, ncol=ncol(x), nrow=ngap), nx[[i]])
  do.call(rbind, nx)
}
