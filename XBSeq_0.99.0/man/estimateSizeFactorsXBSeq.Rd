\name{estimateSizeFactorsXBSeq}
\docType{methods}
\alias{estimateSizeFactorsXBSeq}
\alias{estimateSizeFactorsXBSeq,XBSeqDataSet-method}
\title{Estimate size factors for a XBSeqDataSet object}
\description{
   Same method is adopted from DESeq to estimate the size factors for a XBSeqDataSet
}
\usage{
\S4method{estimateSizeFactorsXBSeq}{XBSeqDataSet}( object, locfunc=median, ... )
}
\arguments{
  \item{object}{
     a XBSeqDataSet
}
  \item{locfunc}{
      a function to compute a location for a sample. 
}   
   \item{...}{extra arguments are ignored}
}
\value{
   The \code{XBSeqDataSet} cds, with the slots \code{sizefactors} updated.
}
\author{
   Yuanhang Liu
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\seealso{
\code{\link{estimateSizeFactorsForMatrixXBSeq}}
}
\examples{
   cds <- makeExampleXBSeqDataSet()
   cds <- estimateSizeFactorsXBSeq( cds )
   sizeFactors( cds )
}