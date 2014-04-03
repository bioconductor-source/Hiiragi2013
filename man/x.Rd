\name{x}
\alias{x}
\docType{data}
\title{Normalized Microarray Data}
\description{
An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object containing the RMA normalized dataset in the assayData and annotation in the phenoData. 
}
\usage{data(x)}
\format{\preformatted{
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
  .. .. ..@ preprocessing    :List of 2
  .. .. .. ..$ filenames  : chr [1:101] "/tmp/RtmpI5T6yy/1_C32_IN.CEL" "/tmp/RtmpI5T6yy/2_C32_IN.CEL" "/tmp/RtmpI5T6yy/3_C32_IN.CEL" "/tmp/RtmpI5T6yy/4_C32_IN.CEL" ...
  .. .. .. ..$ affyversion: chr NA
  .. .. ..@ other            :List of 1
  .. .. .. ..$ : chr ""
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ assayData        :<environment: 0x167aa270> 
  ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':  5 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:5] NA NA NA NA ...
  .. .. ..@ data             :'data.frame':	101 obs. of  5 variables:
  .. .. .. ..$ File.name            : chr [1:101] "1_C32_IN" "2_C32_IN" "3_C32_IN" "4_C32_IN" ...
  .. .. .. ..$ Embryonic.day        : Factor w/ 3 levels "E3.25","E3.5",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ Total.number.of.cells: Factor w/ 11 levels "32","33","34",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ lineage              : chr [1:101] "" "" "" "" ...
  .. .. .. ..$ genotype             : Factor w/ 2 levels "FGF4-KO","WT": 2 2 2 2 2 2 2 2 2 2 ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	3 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr [1:3] NA NA NA
  .. .. ..@ data             :'data.frame':	45101 obs. of  3 variables:
  .. .. .. ..$ symbol  : chr [1:45101] "Copg1" "Atp6v0d1" "Golga7" "Psph" ...
  .. .. .. ..$ genename: chr [1:45101] "coatomer protein complex, subunit gamma 1" "ATPase, H+ transporting, lysosomal V0 subunit D1" "golgi autoantigen, golgin subfamily a, 7" "phosphoserine phosphatase" ...
  .. .. .. ..$ ensembl : chr [1:45101] "ENSMUSG00000030058" "ENSMUSG00000013160" "ENSMUSG00000015341" "ENSMUSG00000029446" ...
  .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ annotation       : chr "mouse4302"
  ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. ..$ labelDescription: chr NA
  .. .. ..@ data             :'data.frame':	101 obs. of  1 variable:
  .. .. .. ..$ ScanDate: chr [1:101] "2011-03-16T04:33:05Z" "2011-03-16T04:42:32Z" "2011-03-16T04:51:54Z" "2011-03-16T05:19:47Z" ...
  .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. ..$ : int [1:3] 1 1 0
  ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slots
  .. .. ..@ .Data:List of 4
  .. .. .. ..$ : int [1:3] 3 0 2
  .. .. .. ..$ : int [1:3] 2 22 0
  .. .. .. ..$ : int [1:3] 1 3 0
  .. .. .. ..$ : int [1:3] 1 0 0
}}
\references{Ohnishi et al., 2014}
\seealso{\code{\link{a}}, \code{\link{xq}}, \code{\link{xql}}}
\examples{
data(x)
x
pData(x)
}
\keyword{datasets}
