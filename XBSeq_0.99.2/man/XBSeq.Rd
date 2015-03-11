\name{XBSeq}
\alias{XBSeq}
\title{
   Express function to carry out XBSeq analysis
}
\description{
   A wrapper function to carry out XBSeq analysis procedure
}
\usage{
XBSeq(observe, background, conditions, method = "pooled", sharingMode = "maximum", fitType = "local", pvals_only = FALSE)
}
\arguments{
  \item{observe}{
   A data matrix contains the observed signal 
}
  \item{background}{
   A data matrix contains the background noise
}
  \item{conditions}{
   A factor to specify the experimental design 
}
  \item{method}{
   Method used to estimate SCV
}
  \item{sharingMode}{
   Mode of sharing of information
}
  \item{fitType}{
   Option to fit mean-SCV relation
}
  \item{pvals_only}{
   Logical; Specify whether to extract pvalues only
}
}
\value{
\item{id}{
   rownames of XBSeqDataSet
}
\item{baseMean}{
   The basemean for all genes
}
\item{baseMeanA}{
   The basemean for condition 'A'
}
\item{baseMeanB}{
   The basemean for condition 'B'
}
\item{foldChange}{
   The fold change compare condition 'B' to 'A'
}
\item{log2FoldChange}{
   The log2 fold change
}
\item{pval}{
   The p value for all genes
}
\item{padj}{
   The adjusted p value for all genes
}
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{estimateRealcount}},  \code{\link{XBSeqDataSet}}, \code{\link{estimateSCV}}, \code{\link{XBSeqTest}}
}
\examples{
   set.seed(1990)
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   colnames(observe_signal) <- paste('S', 1:4, sep='')
   rownames(observe_signal) <- paste('G', 1:5000, sep='')
   colnames(background_noise) <- paste('S', 1:4, sep='')
   rownames(background_noise) <- paste('G', 1:5000, sep='')
   Stats <- XBSeq(observe_signal,background_noise, factor(c('C1','C1','C2','C2') ))
}