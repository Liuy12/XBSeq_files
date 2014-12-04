\name{newXBSeqDataSet}
\alias{newXBSeqDataSet}
\title{
   Create a XBSeqDataSet object
}
\description{
   This function creates a XBSeqDataSet object from a matrix or data frame of count data.
}
\usage{
newXBSeqDataSet(countData, conditions, sizeFactors = NULL, phenoData = NULL, featureData = NULL)
}
\arguments{
  \item{countData}{
      A matrix or data frame of count data, each row indicates each gene; each column indicates each sample
}
  \item{conditions}{
      A factor of experimental conditions (or treatments, or tissue types, or phenotypes, or the like).
      The length of the factor has to be equal to the number of columns of the countData matrix.
}
  \item{sizeFactors}{
      This argument is deprecated. Do not use it. 
}
  \item{phenoData}{
      You may pass an AnnotatedDataFrame here to describe the columns of the count matrix. Note that the
      package always adds two rows (or creates a new AnnotatedDataFrame with only these two rows in case
      you do not supply one) with names "condition" and "sizeFactor" to store this information.
}
  \item{featureData}{
      You may pass an AnnotatedDataFrame here to describe the rows of the count matrix. The package will just
      pass through this information without using it. Note that further columns will be added to feature data
      later, when estimating dispersions.
}
}
\details{
   See also \code{\link{CountDataSet-class}} and the documentation of \code{eSet} (package Biobase) for the meaning of the other slots, which XBSeqDataSet
   inherits from eSet (but which the present package does not use).
}
\value{
   an object of class XBSeqDataSet
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
countsTable <- counts( makeExampleXBSeqDataSet() )
cds <- newXBSeqDataSet( countsTable, c( "A", "A", "A", "B", "B" ) )
}
