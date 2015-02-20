\name{estimateRealcount}
\alias{estimateRealcount}
\alias{insertRow}
\title{
   Preliminary step to estimate the underneath true signal based on observed signal
   and background noise
}
\description{
   Based on the observed signal as well as the background noise, estimate the true signal 
   for each gene. 
}
\usage{
   estimateRealcount(observe, background)
}
\arguments{
  \item{observe}{
   A data.frame which contains the observed signal information; Each row indicates each gene
   and each column indicates each sample.
}
  \item{background}{
   A data.frame which contains the background noise information; The rownames should be the
   same type of annotation as \code{observe}, such as gene symbols, refseqs, etc. 
}
}
\details{
   The observed signal can be achieved by using HTSeq to count the reads map to exonic regions.
   The background noise can be extracted by using HTSeq the second time to count the reads map
   to non-exonic regions, the regions we defined by excluding potential functional elements. 
   The the underneath true signal is estimated by the simple subtraction of observed signal and
   background noise. The true signal of genes with background noise larger than observed signal
   will be assigned as 0.  
}
\value{
   A data.frame contains the estimated true signal for each gene with the same length as observed
   signal.
}
\author{
   Yuanhang Liu
}
\examples{
   set.seed(1990)
   observe_signal <- matrix(rnbinom(20000,100,0.5),nrow=5000,ncol=4)
   background_noise <- matrix(rpois(20000,10),nrow=5000,ncol=4)
   colnames(observe_signal) <- paste('S', 1:4, sep='')
   rownames(observe_signal) <- paste('G', 1:5000, sep='')
   colnames(background_noise) <- paste('S', 1:4, sep='')
   rownames(background_noise) <- paste('G', 1:5000, sep='')   
   Signal <- estimateRealcount(observe_signal,background_noise)
}
