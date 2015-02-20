\name{getsignalVars}
\alias{getsignalVars}
\title{
   Estimate variance of the signal based on variance summation law
}
\description{
   Based on variance of observed signal as well as true signal, estimate the variance of the true signal
}
\usage{
getsignalVars(counts, bgcounts)
}
\arguments{
  \item{counts}{
   count data for observed signal
}
  \item{bgcounts}{
   count data for background noise
}
}
\details{
   By assuming that the true signal and background noise are independent, the variance of the underneath signal (\eqn{\sigma_s^2})
   can be estimated by applying variance summation law:
   \deqn{\sigma_s^2 = \sigma_x^2 + \sigma_b^2 - 2{\rho}{\sigma_x}{\sigma_b}}
   where \eqn{\sigma_x^2} and \eqn{\sigma_b^2} are variance for observed signal and background noise respectively. 
}
\value{
   The estimated variance for true signal
}
\author{
   Yuanhang Liu
}
\seealso{
\code{\link{estimateSCV}}
}
\examples{
   set.seed(1990)
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   colnames(observe_signal) <- paste('S', 1:4, sep='')
   rownames(observe_signal) <- paste('G', 1:5000, sep='')
   colnames(background_noise) <- paste('S', 1:4, sep='')
   rownames(background_noise) <- paste('G', 1:5000, sep='')
   data_var <- getsignalVars(observe_signal, background_noise)  
}