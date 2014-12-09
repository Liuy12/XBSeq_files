\name{plotSCVEsts}
\alias{plotSCVEsts}
\title{
   Plot estimated SCV 
}
\description{
   Plot estimated SCV based on ggplot2
}
\usage{
plotSCVEsts(XB, name = NULL, ymin, linecol = "red3", xlab = "mean of normalized counts", ylab = "SCV")
}
\arguments{
  \item{XB}{
   A XBSeqDataSet object
}
  \item{name}{
   The name of the fit infromation. Only specify this if you choose \code{method="per-condition"}
}
  \item{ymin}{
   The limit of y axis
}
  \item{linecol}{
   The linecolour of the SCV-mean trend
}
  \item{xlab}{
   The lable of x axis
}
  \item{ylab}{
   The lable of y axis
}
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
   Signal <- estimateRealcount(observe_signal, background_noise)
   XB <- newXBSeqDataSet(Signal,factor(c('C1','C1','C2','C2')))
   XB <- estimateSizeFactorsXBSeq( XB )
   XB <-estimateSCV( XB, observe_signal, background_noise)
   plotSCVEsts(XB)
}