\name{dispTable}
\alias{dispTable}
\alias{dispTable,XBSeqDataSet-method}
\alias{dispTable<-}
\alias{dispTable<-,XBSeqDataSet-method}
\title{
    Access the dispersion information for a XBSeqDataSet object
}
\description{
   A method adopted from DESeq to examine the dispersion information for a XBSeqDataSet object
}
\usage{
dispTable(object, ...)
}
\arguments{
  \item{object}{
      a XBSeqDataSet
}
  \item{...}{
      further argumnts are ignored
}  
}
\value{
   Dispersion information for a XBSeqDataSet object
}
\author{
   Simon Anders
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\seealso{
\code{\link{estimateSCV}}, \code{\link{dispEst}}, \code{\link{XBSeqDataSet}}
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
   dispTable( XB )
}   
