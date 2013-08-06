myHeatmap2 <-
function(x, 
  rowGroups = factor(rep(1, nrow(x))),
  colGroups = factor(rep(1, ncol(x))),
  keeprownames  = TRUE,
  colors = colorRampPalette(brewer.pal(9, "Blues")[-1])(100),
  ...) {

  x = ordermat(x,    rowGroups, ...)
  x = ordermat(t(x), colGroups, ...)
  
  if(!keeprownames)
    colnames(x) = NULL

  print(levelplot(x,
          aspect = "fill",
          xlab="", ylab="",
          scales = list(x = list(rot = 90), raster=TRUE),
          col.regions = colors, 
          colorkey = list(space = "left", height=0.15, useRaster=TRUE)))
}
