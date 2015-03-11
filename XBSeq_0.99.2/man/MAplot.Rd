\name{MAplot}
\alias{MAplot}
\title{
   Generate maplot after differential expression test
}
\description{
   Generate maplot after differential expression test based on ggplot2
}
\usage{
MAplot(stats, ylim, padj = TRUE, pcuff = 0.1, lfccuff = 1, linecol = "red3", xlab = "mean of normalized counts", ylab = expression(log[2] ~ fold ~ change), shape) 
}
\arguments{
  \item{stats}{
   The output of \code{XBSeqTest}
}
  \item{ylim}{
   Range of limit for y axis
}
  \item{padj}{
   Whether to use adjusted p value or not
}
  \item{pcuff}{
   Threshold for pvalue
}
  \item{lfccuff}{
   Log fold change cutoff
}
  \item{linecol}{
   Colour of horizontal line 
}
  \item{xlab}{
   Lable for x axis
}
  \item{ylab}{
   Lable for y axis
}
  \item{shape}{
   The shape of the points used
}  
}
\details{
   Generate classic MAplot for DE analysis using ggplot2, where A and M are from slot \code{baseMean} and slot \code{log2FoldChange} of 
   the test statistics aftering calling \code{XBSeqTest}. The ggplot2 package generally generate figures of better quality as well as give 
   user better control of the plotting system compared with the base plotting system. 
}
\value{
   MAplot of test statistics
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{\code{\link{XBSeqTest}}}
\examples{
   set.seed(1990)
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   colnames(observe_signal) <- paste('S', 1:4, sep='')
   rownames(observe_signal) <- paste('G', 1:5000, sep='')
   colnames(background_noise) <- paste('S', 1:4, sep='')
   rownames(background_noise) <- paste('G', 1:5000, sep='')
   stats <- XBSeq (observe_signal,background_noise,factor(c('C1','C1','C2','C2')) )
   MAplot(stats)
}