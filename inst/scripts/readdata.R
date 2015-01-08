library("affy")
library("ArrayExpress")
library("arrayQualityMetrics")
library("mouse4302.db")
library("RColorBrewer")

CELdir    = tempdir()
CELfiles  = getAE("E-MTAB-1681", path = CELdir, type = "raw")$rawFiles

## -------------------------------------------------------------------------------------
## Read array metadata table and fill empty cells in the columns Embryonic.day
## and Total.number.of.cells by the values implied ## by the non-empty cells above
## -------------------------------------------------------------------------------------

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
  
  row.names(x) = paste(x$File.name, "CEL", sep = ".")
  return(x)
}

## -------------------------------- Script starts here ---------------------------------

pdata = readCSVtable(system.file("scripts", "annotation.csv", package = "Hiiragi2013"))

pdata$genotype = as.factor(ifelse(grepl("_KO$", pdata$File.name), "FGF4-KO", "WT"))

## -------------------------------------------------------------------------------------
## Read the CEL files
## -------------------------------------------------------------------------------------

fileNames = row.names(pdata)
fileExists = (fileNames %in% CELfiles)
stopifnot(all(fileExists))

a = ReadAffy(filenames = fileNames, celfile.path = CELdir, phenoData = pdata,
             verbose = TRUE)

pData(a)$ScanDate = factor(as.Date(sub( "10/16/09", "2010-09-16",
    sapply(strsplit( protocolData(a)$ScanDate, split = "[T ]" ), `[`, 1) )))

save(a, file="a.rda", compress="xz")

## -------------------------------------------------------------------------------------
## Normalize with RMA
## -------------------------------------------------------------------------------------

x = rma(a)

## Create columns
## fData(x)$symbol:   gene symbols where available, Affy feature ID otherwise
## fData(x)$genename: a more verbose gene description

annotateGene = function(db, what, missing) {
  tab = toTable(db[ featureNames(x) ])
  mt = match( featureNames(x), tab$probe_id)
  ifelse(is.na(mt), missing, tab[[what]][mt])
}
fData(x)$symbol   = annotateGene(mouse4302SYMBOL,   "symbol", missing = featureNames(x))
fData(x)$genename = annotateGene(mouse4302GENENAME, "gene_name", missing = "")
fData(x)$ensembl  = annotateGene(mouse4302ENSEMBL,  "ensembl_id", missing = "")

## -------------------------------------------------------------------------------------
## Grouping of samples 
## -------------------------------------------------------------------------------------

## We define a grouping of the samples and an associated colour map, which will
## be used in the plots throughout this report.

groups = with(pData(x), list(
  `E3.25`           = which(genotype=="WT" & Embryonic.day=="E3.25"),
  `E3.5 (EPI)`      = which(genotype=="WT" & Embryonic.day=="E3.5" & lineage=="EPI"),
  `E4.5 (EPI)`      = which(genotype=="WT" & Embryonic.day=="E4.5" & lineage=="EPI"),
  `E3.5 (PE)`       = which(genotype=="WT" & Embryonic.day=="E3.5" & lineage=="PE"),
  `E4.5 (PE)`       = which(genotype=="WT" & Embryonic.day=="E4.5" & lineage=="PE"),
  `E3.25 (FGF4-KO)` = which(genotype=="FGF4-KO" & Embryonic.day=="E3.25"),
  `E3.5 (FGF4-KO)`  = which(genotype=="FGF4-KO" & Embryonic.day=="E3.5"),
  `E4.5 (FGF4-KO)`  = which(genotype=="FGF4-KO" & Embryonic.day=="E4.5")))

sampleColourMap = character(length(groups))
names(sampleColourMap) = names(groups)
sampleColourMap[c("E3.5 (EPI)", "E4.5 (EPI)")]         = brewer.pal(10, "Paired")[1:2]
sampleColourMap[c("E3.5 (PE)", "E4.5 (PE)")]           = brewer.pal(10, "Paired")[3:4]
sampleColourMap[c("E3.25 (FGF4-KO)")]                  = brewer.pal(10, "Paired")[c(7)]
sampleColourMap[c("E3.5 (FGF4-KO)", "E4.5 (FGF4-KO)")] = brewer.pal(10, "Paired")[c(8,6)]
sampleColourMap[c("E3.25")]                            = brewer.pal(12, "Paired")[c(9)]
stopifnot(!any(sampleColourMap==""))

## The following assertions aim to make sure that each sample was assigned to
## exactly one group.
stopifnot(!any(duplicated(unlist(groups))),
            setequal(unlist(groups), seq_len(ncol(x))),
            setequal(names(sampleColourMap), names(groups)))

## Next, assign a colour and a name to each sample, which will be used in the
## subsequent plots. For sample names, use the group name and the array numeric
## index.

sampleNames = sampleGroup = rep(NA_character_, ncol(x))
for(i in seq(along = groups)) {
  idx = groups[[i]]
  sampleGroup[idx] = names(groups)[i]
  sampleNames[idx] = paste(idx, names(groups)[i])
}
pData(x)$sampleGroup  = sampleGroup
pData(x)$sampleColour = sampleColourMap[sampleGroup]
sampleNames(x) = sampleNames

save(x, file="x.rda", compress="xz")
