MDSplot <-
function(x, mask,
     flip=integer(0), rotation = 0,
     cex=2, col=x$sampleColour,
     panellabel, pointlabel) {
  
  mat = exprs(x)
  distances = dist(t(mat))
  mds = isoMDS(distances)$points
  for(i in flip)
    mds[,i] = -mds[,i]
  mds = mds %*% rotmat(rotation)
  
  par(mai=c(0.1, 0.1, 0.1, 0.1))
  xlim = range(mds[,1])
  ylim = range(mds[,2])
  if(missing(mask))
      mask = rep(TRUE, ncol(x))
  theColours = paste(col[mask], "a0", sep="")
  plot(mds[mask,], pch = 16,
    cex=cex, xlab = "", ylab = "", xaxt = "n", yaxt="n",
    asp = 1, xlim = xlim, ylim = ylim,
    col = if(missing(pointlabel)) theColours else par("bg"))
  if(missing(pointlabel)) {
    points(mds[mask,], pch = 1,cex=cex, col=grey(0.3)) 
  } else {
    text(mds[mask,], labels=pointlabel,
         adj=c(0.5, 0.5), cex=cex*2/3, col=theColours)
  }
  if(!missing(panellabel))
      text(par("usr")[1], par("usr")[4], panellabel,
       adj=c(-1, 1.5), cex=cex*2/3, col="#808080")
}
