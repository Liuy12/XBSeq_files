\name{getCountParams}
\alias{getCountParams}
\title{
   Extract the mean and variance for a count dataset
}
\description{
   Internal function to extract the mean and variance for each gene based on the normalized counts
}
\usage{
getCountParams(counts, sizeFactors)
}
\arguments{
  \item{counts}{
   one integer count data
}
  \item{sizeFactors}{
   sizeFactors calculated by \code{estimateSizeFactorsXBSeq}
}
}
\value{
   a data.frame that contains the mean and variance calculated for each gene
}
\author{
   Yuanhang Liu
}
\examples{
   cds <- makeExampleXBSeqDataSet()
   cds <- estimateSizeFactorsXBSeq( cds )
   getCountParams(counts(cds),sizeFactors(cds))
}