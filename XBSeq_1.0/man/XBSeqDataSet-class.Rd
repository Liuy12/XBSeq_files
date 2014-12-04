\name{XBSeqDataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{XBSeqDataSet-class}
\alias{conditions,XBSeqDataSet-method}
\alias{conditions<-,XBSeqDataSet-method}
\alias{counts,XBSeqDataSet-method}
\alias{counts<-,XBSeqDataSet,matrix-method}
\alias{dispTable,XBSeqDataSet-method}
\alias{dispTable<-,XBSeqDataSet-method}
\alias{estimateSCV,XBSeqDataSet-method}
\alias{estimateSizeFactorsXBSeq,XBSeqDataSet-method}
\alias{sizeFactors,XBSeqDataSet-method}
\alias{sizeFactors<-,XBSeqDataSet,numeric-method}

\title{Class \code{"XBSeqDataSet"}}
\description{
   Package for use in XBSeq
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("XBSeqDataSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fitInfo}:}{Object of class \code{"environment"} ~~ }
    \item{\code{dispTable}:}{Object of class \code{"character"} ~~ }
    \item{\code{multivariateConditions}:}{Object of class \code{"logical"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{conditions}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{conditions<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{counts}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{counts<-}{\code{signature(object = "XBSeqDataSet", value = "matrix")}: ... }
    \item{dispTable}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispTable<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateSCV}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateSizeFactorsXBSeq}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{sizeFactors}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{sizeFactors<-}{\code{signature(object = "XBSeqDataSet", value = "numeric")}: ... }
	 }
}
\author{
   Yuanhang Liu
}