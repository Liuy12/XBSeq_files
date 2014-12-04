\name{getSCV}
\alias{getSCV}
\title{
   Extract estimation of squared coefficient of variation 
}
\description{
   A internal function called by \code{estimateSCV}. There is no need to call this function directly
}
\usage{
getSCV(means, variances, sizeFactors, fitType = c("parametric", "local"), locfit_extra_args = list(), lp_extra_args = list(), adjustForBias = TRUE)
}
\arguments{
  \item{means}{
   Mean statistics for each gene
}
  \item{variances}{
   Variance statistics for each gene
}
  \item{sizeFactors}{
   sizeFactors for a XBSeqDataSet object
}
  \item{fitType}{
   The method that will be used to fit mean-SCV relation, can either be 'parametric', or 'local'
}
  \item{locfit_extra_args}{
   Futher arguments supplied for "locfit"
}
  \item{lp_extra_args}{
   Futher arguments supplied for \code{lp} function in \code{locfit} package
}
  \item{adjustForBias}{
   Logical; set whether to carry out bias correction procedure for SCV
}
}
\value{
   a fit object based on means and SCV
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
\code{\link{estimateSCV}}
}
\examples{
   cds <- makeExampleXBSeqDataSet()
   data <- getCountParams( counts(cds), sizeFactors(cds), conditions(cds) )
   SCVf <- getSCV( data$baseMean, data$baseVar, sizeFactors(cds))
}