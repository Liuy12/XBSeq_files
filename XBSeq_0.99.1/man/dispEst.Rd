\name{dispEst}
\docType{methods}
\alias{dispEst}
\alias{dispEst,XBSeqDataSet-method}
\alias{dispEst<-}
\alias{dispEst<-,XBSeqDataSet-method}
\title{
   Function to access the dispersion estimation for each gene 
}
\description{
   The dispersion estimated for each gene are stored as a data.frame after user called \code{\link{estimateSCV}}
}
\usage{
dispEst(object, varname = NA)
dispEst(object, varname = NA) <- value
}
\arguments{
  \item{object}{
   XBSeqDataSet object
}
  \item{varname}{
   variable name of dispersion estimates
}
  \item{value}{
   The dispersion estimates for each gene 
}
}
\value{
   A data.frame which contains the dispersion estimates for each gene
}
\author{
   Yuanhang Liu
}
\seealso{
\code{\link{estimateSCV}}, \code{\link{dispTable}}, \code{\link{XBSeqDataSet}}
}
\examples{
   set.seed(1990)
   conditions <- factor(c('C1','C1','C2','C2'))
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   colnames(observe_signal) <- paste('S', 1:4, sep='')
   rownames(observe_signal) <- paste('G', 1:5000, sep='')
   colnames(background_noise) <- paste('S', 1:4, sep='')
   rownames(background_noise) <- paste('G', 1:5000, sep='')
   Signal <- estimateRealcount(observe_signal, background_noise)
   XB <- XBSeqDataSet(Signal,conditions)
   XB <- estimateSizeFactors( XB )
   XB <- estimateSCV( XB, observe_signal, background_noise,fitType='local')
   dispEst( XB )
}  
