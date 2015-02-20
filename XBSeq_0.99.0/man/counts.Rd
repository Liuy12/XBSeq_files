\name{counts}
\docType{methods}
\alias{counts}
\alias{counts,XBSeqDataSet-method}
\alias{counts<-,XBSeqDataSet,matrix-method}
\title{
   Access the counts table in XBSeqDataSet object.
}
\description{
   The counts slot holds the count data as a matrix of non-negative integer
   count values, similar method adopted from DESeq
}
\usage{
\S4method{counts}{XBSeqDataSet}(object, normalized=FALSE)
\S4method{counts}{XBSeqDataSet,matrix}(object) <- value
}
\arguments{
  \item{object}{a \code{XBSeqDataSet} object.}
  \item{normalized}{logical indicating whether or not to divide the
    counts by the size factors before returning.}
  \item{value}{the integer count matrix for a XBSeqDataSet.}  
}
\value{
   Either the count matrix or the normalized matrix if the argument normalized
   is set to TRUE 
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
   head( counts( cds ) )
}   
   
