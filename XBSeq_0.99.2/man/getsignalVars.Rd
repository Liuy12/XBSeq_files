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
   A data frame or matrix which contains the expression level information for an experiment. Rows represent genes and 
   Columns represent samples. Please refer details for more information. 
}
  \item{bgcounts}{
   A data frame or matrix which contains the background intensities information for an experiment. Rows represent genes and 
   Columns represent samples. Please refer details for more information. 
}
}
\details{
   Observed signal are the reads mapped to the exonic regions which can be obtained by applying htseq-count procedure with GTF files of 
   exonic regions. Background noise are the reads mapped to the non-exonic regions which can be obtained by applying htseq-count procedure 
   with GTF files of non-exonic regions we defined by certain criteria. Details regarding how to carry out the htseq procedure for observed
   signal as well as background noise can be found in the vignette of XBSeq. One example dataset is provided in \code{\link{ExampleData}}. 
   
   By assuming that the true signal and background noise are independent, the variance of the underneath signal (\eqn{\sigma_s^2})
   can be estimated by applying variance summation law:
   \deqn{\sigma_s^2 = \sigma_x^2 + \sigma_b^2 - 2{\rho}{\sigma_x}{\sigma_b}}
   where \eqn{\sigma_x^2} and \eqn{\sigma_b^2} are variance for observed signal and background noise respectively. 
}
\value{
   A data frame or matrix with the same number of rows as \code{counts}. Rows represent the estimated variance of true signal for each gene. 
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{estimateSCV}}
}
\examples{
   set.seed(1990)
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   colnames(observe_signal) <- paste('S', 1:4, sep='')
   rownames(observe_signal) <- paste('G', 1:5000, sep='')
   colnames(background_noise) <- paste('S', 1:4, sep='')
   rownames(background_noise) <- paste('G', 1:5000, sep='')
   data_var <- getsignalVars(observe_signal, background_noise)  
}