\name{fitInfo}
\alias{fitInfo}
\title{
    Accessor function for the fitInfo objects in a XBSeqDataSet
}
\description{
   Same method is adopted from DESeq to access the fit information from a XBSeqDataSet
}
\usage{
fitInfo( XB, name=NULL )
}
\arguments{
  \item{XB}{
      a XBSeqDataSet
}
  \item{name}{
      if \code{estimateSCV} was called with \code{method="per-condition"}
      a name hasd to specified. Try \code{ls(XB@fitInfo}.
}
}
\seealso{\code{\link{estimateSCV}} }
\author{
   Yuanhang Liu
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\examples{
  set.seed(1990)
   conditions <- factor(c('C1','C1','C2','C2'))
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   Signal <- estimateRealcount(observe_signal, background_noise)
   XB <- newXBSeqDataSet(Signal,conditions)
   XB <- estimateSizeFactorsXBSeq( XB )
   XB <- estimateSCV( XB, observe_signal, background_noise,fitType='local')
   str( fitInfo( XB ) )
}   

