\name{adjustScv}
\alias{adjustScv}
\title{
   Adjust bias for estimating coefficient of variation (scv)
}
\description{
   The same method from DESeq is adopted to adjust the bias for SCV. We carried out 
   bias correction procedure at estimated signal level rather than observed signal level as in DESeq  
}
\usage{
adjustScv(scv, nsamples)
}
\arguments{
  \item{scv}{
   Raw squared coefficient of variation (SCV) for estimated true signal.
}
  \item{nsamples}{
   The number of replicates in each condition
}
}
\value{
   An unbiased estimate for SCV
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
# This example is adopted from DESeq
   true_mean <- 100
   true_scv <- .1
   nsamples <- 3
   res <- replicate( 1, {
      mySample <- rnbinom( nsamples, mu=true_mean, size=1/true_scv )
      mu_est <- mean( mySample )
      raw_var_est <- var( mySample ) - mean( mySample )
      raw_scv_est <- raw_var_est / mu_est^2
      unbiased_raw_scv_est <- adjustScv( raw_scv_est, 4 )
      c( raw_scv_est = raw_scv_est, unbiased_raw_scv_est = unbiased_raw_scv_est ) } )
   rowMeans( res )
}