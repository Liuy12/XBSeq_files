\name{makeExampleXBSeqDataSet}
\alias{makeExampleXBSeqDataSet}
\title{
   Generate an example XBSeqDataSet object
}
\description{
   The same function is adopted from DESeq. This function returns an example XBSeqDataSet. 
   It is used for the examples in the package help pages.
}
\usage{
makeExampleXBSeqDataSet()
}
\value{
   a XBSeqDataSet that has been constructed as follows: First, true base mean
   values for 10,000 genes are drawn from an exponential distribution with rate
   1/250. Then, certain genes are declared (with probability 0.3 per gene) as
   truly differentially expressed (tDE). For these genes, the true base mean
   is split into two values, one for condition "A" and one for condition "B",
   such that the log2 fold change from "A" to "B" follows a zero-centred
   normal distribution with standard deviation 2. Then, counts are drawn for
   each gene for 5 samples, the first three corresponding to condition "A" and
   the remaining two for condition "B". The counts are drawn from a negative
   binomial with the specified mean, multiplied by the size factor for the sample,
   with a constant raw SCV (dispersion) of 0.2 (i.e., a 'size' parameter of 1/0.2). The true
   size factors are fixed to c( 1., 1.3, .7, .9, 1.6 ). 
}
\author{
   Yuanhang Liu
}
\references{
   Simon Anders, Wolfgang Huber: Differential expression
   analysis for sequence count data. Genome Biology 11 (2010)
   R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
}
\examples{
   data <- makeExampleXBSeqDataSet
}
