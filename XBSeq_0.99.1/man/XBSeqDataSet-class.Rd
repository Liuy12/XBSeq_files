\name{XBSeqDataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{XBSeqDataSet-class}
\alias{XBSeqDataSet}
\alias{XBSeqDataSetFromHTSeqCount}
\alias{makeExampleXBSeqDataSet}

\title{Class \code{"XBSeqDataSet"}}
\description{
   Class to be used in XBSeq package. \code{"XBSeqDataSet"} is inherited from \code{"DESeqDataSet"}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("XBSeqDataSet", XB, fitInfo, dispTable, conditions, dispEst, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fitInfo}:}{Object of class \code{"environment"} }
    \item{\code{dispTable}:}{Object of class \code{"character"} }
    \item{\code{conditions}:}{Object of class \code{"character"} }
    \item{\code{dispEst}:}{Object of class \code{"data.frame"} }
    \item{\code{assays}:}{A list or SimpleList of matrix elements, or a matrix as used in class \code{"SummarizedExperiment"}}
    \item{\code{colData}:}{An optional DataFrame describing the samples as used in class \code{"SummarizedExperiment"} }
    \item{\code{design}:}{a formula which specifies the design of the experiment, taking the form \code{formula(~ x + y + z)} }
    \item{\code{dispersionFunction}:}{The dispersion function as used in class \code{"DESeqDataSet"}  }
    \item{\code{exptData}:}{An optional \code{"SimpleList"} of arbitrary content describing the overall experiment as used in class \code{"SummarizedExperiment"} }
    \item{\code{rowData}:}{A GRanges or GRangesList instance describing the ranges of interest as used in class \code{"SummarizedExperiment"} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{DESeqDataSet}"}, directly
}
\section{Methods}{
  \describe{
    \item{conditions}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{conditions<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{counts}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{counts<-}{\code{signature(object = "XBSeqDataSet", value = "matrix")}: ... }
    \item{dispTable}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispTable<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispEst}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispEst<-}{\code{signature(object = "XBSeqDataSet")}: ... }
	\item{fitInfo}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{fitInfo<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateSCV}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateSizeFactorsXBSeq}{\code{signature(object = "XBSeqDataSet")}: ... }
	 }
}
\author{
   Yuanhang Liu
}