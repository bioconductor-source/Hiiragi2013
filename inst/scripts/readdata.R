library("affy")
library("ArrayExpress")
library("arrayQualityMetrics")
library("mouse4302.db")

CELdir    = tempdir()
CELfiles  = getAE("E-MTAB-1681", path = CELdir, type = "raw")$rawFiles

## ---------------------------------------------------------------------------------
## Read array metadata table and fill empty cells in the columns Embryonic.day
## and Total.number.of.cells by the values implied ## by the non-empty cells above
## ---------------------------------------------------------------------------------

fillColumn = function(x, empty){
  wh  = which(!empty(x))
  len = length(wh)
  wh  = c(wh, length(x)+1)
  for(i in seq_len(len))
    x[ wh[i]:(wh[i+1]-1) ] = x[wh[i]]
  return(x)
}

readCSVtable = function(name) {
  x = read.csv(name, stringsAsFactors = FALSE, colClasses = "character")
  x$Embryonic.day = factor(fillColumn(x$Embryonic.day, empty = function(x) x==""))
  
  wh = which(colnames(x)=="Total.number.of.cells")
  if(length(wh)==1) {
    x[[wh]] = fillColumn(x$Total.number.of.cells, empty = function(x) x=="" & !is.na(x))
    x[[wh]] = as.integer(x[[wh]])
  } else {
    x$"Total.number.of.cells" = rep(NA, nrow(x))
  }
  x$Total.number.of.cells = addNA(as.factor(x$Total.number.of.cells))
  
  wh = which(colnames(x) %in% c("X.EPI...PE.", "Type"))
  if(length(wh)==1) {
    colnames(x)[wh] = "lineage"
  } else {
    x$lineage = rep(NA, nrow(x))
  }
  
  return(x)
}

## ------------------------------ Script starts here -------------------------------

pdata = readCSVtable(system.file("scripts", "annotation.csv", package = "Hiiragi2013"))

pdata$genotype = as.factor(ifelse(grepl("_KO$", pdata$File.name), "FGF4-KO", "WT"))

## ---------------------------------------------------------------------------------
## Read the CEL files
## ---------------------------------------------------------------------------------

fileNames = paste(pdata$File.name, "CEL", sep = ".")
fileExists = (fileNames %in% CELfiles)
stopifnot(all(fileExists))

a = ReadAffy(filenames = fileNames, celfile.path = CELdir, phenoData = pdata,
             verbose = TRUE)

pData(a)$ScanDate = factor(as.Date(sub( "10/16/09", "2010-09-16",
    sapply(strsplit( protocolData(a)$ScanDate, split = "[T ]" ), `[`, 1) )))

save(a, file="a.rda")

## ---------------------------------------------------------------------------------
## Normalize with RMA
## ---------------------------------------------------------------------------------

x = rma(a)

## Create columns
## fData(x)$symbol:   gene symbols where available, Affy feature ID otherwise
## fData(x)$genename: a more verbose gene description

annotateGene = function( db, what, missing) {
  tab = toTable(db[ featureNames(x) ])
  mt = match( featureNames(x), tab$probe_id)
  ifelse(is.na(mt), missing, tab[[what]][mt])
}
fData(x)$symbol   = annotateGene(mouse4302SYMBOL,   "symbol", missing = featureNames(x))
fData(x)$genename = annotateGene(mouse4302GENENAME, "gene_name", missing = "")
fData(x)$ensembl  = annotateGene(mouse4302ENSEMBL,  "ensembl_id", missing = "")

save(x, file="x.rda")
