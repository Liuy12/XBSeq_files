\name{conditions}
\docType{methods}
\alias{conditions}
\alias{conditions<-,XBSeqDataSet-method}
\alias{conditions,XBSeqDataSet-method}
\title{
   Accessor functions for the 'conditions' information in a
   XBSeqDataSet object.}
\description{
   Conditions extract the experimental design information similar as used in DESeq. 
}
\usage{
\S4method{conditions}{XBSeqDataSet}(object,...)
\S4method{conditions}{XBSeqDataSet}(object,...) <- value
}
\arguments{
  \item{object}{a XBSeqDataSet}
  \item{value}{experimental design information}
  \item{...}{Further arguments will be ignored}
  }
\value{
   The experimental design information for a XBSeqDataSet object
}
\author{
   Yuanhang Liu
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\examples{
   cds <- makeExampleXBSeqDataSet()
   conditions( cds )
}  