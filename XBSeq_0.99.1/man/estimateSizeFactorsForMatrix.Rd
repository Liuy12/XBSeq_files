\name{estimateSizeFactorsForMatrix}
\alias{estimateSizeFactorsForMatrix}
\title{
   Internal function to estimate the size factors
}
\description{
   The is the same method as used in DESeq and DESeq2 to estimate the size factors
}
\usage{
estimateSizeFactorsForMatrix(counts, locfunc = median)
}
\arguments{
  \item{counts}{
   a matrix of counts information 
}
  \item{locfunc}{
   a function to compute a location for a sample.
}
}
\value{
   Size factors for each sample in the count matrix.
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\author{
   Simon Anders
}
\examples{
   cds <- makeExampleXBSeqDataSet()
   estimateSizeFactorsForMatrix( counts(cds) )
}
