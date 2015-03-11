\name{ExampleData}
\alias{ExampleData}
\alias{Observed}
\alias{Background}
\title{
   Example Datasets used in vignette
}
\description{
   Example datasets used in vignette by carry out HTSeq procedure for exonic mapped reads (Observed) and non-exonic mapped reads (Background).
}
\details{
   In order to use XBSeq for testing DE, we need to run HTSeq twice to measure the reads mapped to exonic regions (observed signal) and 
   non-exonic regions (background noise). Exonic region annotation is obtained from UCSC database. Non-exonic regions are constructed by following
   several criteria: 1) 1) Download refFlat table from UCSC database 2) Several functional elements (mRNA, peudo genes, etc.) of the genome are 
   excluded from the whole genome annotation. 3) Construct the background region for each gene by making the region to have the same sturcture or 
   length as the exonic region of the gene. More details regarding how do we construct the background region annotation file can be found in vignette 
   of XBSeq. 
   
   We carried out HTSeq procedure twice by using a a mouse RNA-seq dataset, which contains 3 replicates of wild type mouse liver tissues (WT) and 3 
   replicates of Myc transgenic mouse liver tissues (MYC). ExampleData contains two datasets. One is called observed. The other one is called background.  
   For the two datasets, rows represent exonic or non-exonic region mapped reads for each gene. Columns represent each sample. The two datasets can be 
   loaded via \code{data(ExampleData)} after loading the XBSeq library. 
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}