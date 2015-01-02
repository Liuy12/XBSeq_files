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
   data_var <- getsignalVars(observe_signal, background_noise)  
}