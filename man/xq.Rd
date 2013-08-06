\name{xq}
\alias{xq}
\docType{data}
\title{
qPCR Gene Expression 
}
\description{
Single-cell gene expression measured by qPCR stored in an \code{ExpressionSet}.
}
\usage{data(xq)}
\format{
  The format is:
Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
  ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
  .. .. ..@ name             : chr ""
  .. .. ..@ lab              : chr ""
  .. .. ..@ contact          : chr ""
  .. .. ..@ title            : chr ""
  .. .. ..@ abstract         : chr ""
  .. .. ..@ url              : chr ""
  .. .. ..@ pubMedIds        : chr ""
  .. .. ..@ samples          : list()
  .. .. ..@ hybridizations   : list()
  .. .. ..@ normControls     : list()
  .. .. ..@ preprocessing    : list()
  .. .. ..@ other            : list()
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ assayData        :<environment: 0x95bc008> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	6 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:6] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	137 obs. of  6 variables:
  .. .. .. ..$ Sample.ID    : chr [1:137] "34c_2" "34c_3" "34c_4" "34c_5" ...
  .. .. .. ..$ Embryonic.day: chr [1:137] "E3.25" "E3.25" "E3.25" "E3.25" ...
  .. .. .. ..$ Cell.type    : chr [1:137] "ICM" "ICM" "ICM" "ICM" ...
  .. .. .. ..$ sampleColours: chr [1:137] "#CAB2D6" "#CAB2D6" "#CAB2D6" "#CAB2D6" ...
  .. .. .. ..$ sampleGroup  : chr [1:137] "E3.25" "E3.25" "E3.25" "E3.25" ...
  .. .. .. ..$ sampleColour : chr [1:137] "#CAB2D6" "#CAB2D6" "#CAB2D6" "#CAB2D6" ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr NA
  .. .. ..@ data             :'data.frame':	38 obs. of  1 variable:
  .. .. .. ..$ symbol: chr [1:38] "Fgf4" "Tom1l1" "Tdgf1" "Cldn4" ...
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ annotation       : chr "single cell qPCR"
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr(0) 
  .. .. ..@ data             :'data.frame':	137 obs. of  0 variables
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 4
  .. .. .. ..$ : int [1:3] 2 15 0
  .. .. .. ..$ : int [1:3] 2 16 0
  .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. ..$ : int [1:3] 1 0 0
}
\references{
Ohnishi et al., 2013
}
\examples{
data(xq)
## maybe str(xq) ; plot(xq) ...
}
\keyword{datasets}
