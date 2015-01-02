\name{sizeFactors}
\docType{methods}
\alias{sizeFactors}
\alias{sizeFactors<-,XBSeqDataSet,numeric-method}
\alias{sizeFactors,XBSeqDataSet-method}
\title{
   Access the sizefactor information from a XBSeqDataSet object
}
\description{
   Access the sizefactors assigned to each samples in a XBSeqDataSet object
}
\usage{
\S4method{sizeFactors}{XBSeqDataSet}(object)
\S4method{sizeFactors}{XBSeqDataSet,numeric}(object) <- value
}
\arguments{
  \item{object}{a \code{DESeqDataSet} object.}
  \item{value}{a numeric vector, one size factor for each column in the count data.}
}
\value{
   sizeFactors of a XBSeqDataSet object
}
\seealso{
\code{\link{estimateSizeFactorsXBSeq}}
}
\examples{
   cds <- makeExampleXBSeqDataSet()
   cds <- estimateSizeFactorsXBSeq( cds )
   sizeFactors(cds)
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\author{
   Yuanhang Liu
}