\name{MAplot}
\alias{MAplot}
\title{
   Generate maplot after differential expression test
}
\description{
   Generate maplot after differential expression test based on ggplot2
}
\usage{
MAplot(stats, ylim, padj = T, pcuff = 0.1, lfccuff = 1, linecol = "red3", xlab = "mean of normalized counts", ylab = expression(log[2] ~ fold ~ change))
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
}
\author{
   Yuanhang Liu
}
\examples{
   set.seed(1990)
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   stats <- XBSeq (observe_signal,background_noise,factor(c('C1','C1','C2','C2')) )
   MAplot(stats)
}