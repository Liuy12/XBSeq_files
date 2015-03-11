\name{XBSeqTest}
\alias{XBSeqTest}
\title{
   XBSeq test for differential expression
}
\description{
   The same method is adopted from DESeq for testing differential expression 
}
\usage{
XBSeqTest(XB, condA, condB, pvals_only = FALSE)
}
\arguments{
  \item{XB}{
   A XBSeqDataSet object
}
  \item{condA}{
   Factor specified for condition A
}
  \item{condB}{
   Factor specified for condition B
}
  \item{pvals_only}{
   Logical;whether or not only extract p values
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
\details{
   Differential expression analysis based on statistical methods proposed for DESeq. Details about the method can be found in DESeq manual page. 
}
   
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\author{
   Yuanhang Liu
}
\seealso{
\code{\link{XBSeq}}, \code{\link{estimateSCV}}
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
   Teststas <- XBSeqTest( XB, levels(conditions)[1L], levels(conditions)[2L])
}