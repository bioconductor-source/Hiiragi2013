getDifferentialExpressedGenes <-
function(x, groups, g1, g2, theta=0.5, FDRcutoff = 0.05) {
   stopifnot(g1 %in% names(groups), g2 %in% names(groups))
   fac = rep(NA_integer_, ncol(x))
   fac[groups[[g1]]] = 1
   fac[groups[[g2]]] = 2
   fac = factor(fac)

   rv  = rowVars(exprs(x)[, !is.na(fac)])
   passfilter = (rv > quantile(rv, probs=theta))

   tt = rowttests(x, fac)

   tt$adjp = rep(NA_real_, nrow(tt))
   tt$adjp[passfilter] = p.adjust(tt$p.value[passfilter])
   ord = order(tt$adjp)
   differentially = ord[ tt$adjp[ord]<=FDRcutoff & !is.na(tt$adjp[ord]) ]
   return(differentially)
 }
