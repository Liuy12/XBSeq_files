\name{dispTable}
\alias{dispTable}
\alias{dispTable,CountDataSet-method}
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
\seealso{\code{\link{estimateSCV}}, \code{\link{XBSeqTest}} }
\author{
   Yuanhang Liu
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\examples{
# Example adopted from DESeq
  set.seed(1990)
   conditions <- factor(c('C1','C1','C2','C2'))
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   Signal <- estimateRealcount(observe_signal, background_noise)
   XB <- newXBSeqDataSet(Signal,conditions)
   XB <- estimateSizeFactorsXBSeq( XB )
   XB <- estimateSCV( XB, observe_signal, background_noise,fitType='local')
   dispTable( XB )
}   
