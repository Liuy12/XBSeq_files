# Introduction 

XBSeq is a novel algorithm for testing RNA-seq differential expression (DE), where a statistical model was established
based on the assumption that observed signals are the convolution of true expression signals and sequencing noises. The
mapped reads in non-exonic regions are considered as sequencing noises, which follows a Poisson distribution. Given
measureable observed and noise signals from RNA-seq data, true expression signals, assuming governed by the negative
binomial distribution, can be delineated and thus the accurate detection of differential expressed genes

# Installation 

Currently we have already finalizing the documentation for XBSeq to be released in bioconductor. Right now the algorithm is implemented in R. To use the algorithm, you can either download the XBSeq_1.0.tar.gz file and install locally or you can install XBSeq by the following steps: 

Firstly, please make sure that the devtools package installed and loaded in your libraries
```r
install.packages("devtools")
library(devtools)
```
Then, install XBSeq directly from github 
```r
install_github("Liuy12/XBSeq/XBSeq_0.99.2")
library('XBSeq')
```
XBSeq depends on Biobase, locfit, pracma, matrixStats, ggplot2 from Bioconductor

# Use XBSeq for testing differential expression 

### HTseq counting

In order to use XBSeq for testing DE, we need to run HTSeq twice to measure the reads mapped to exnoic regions (observed signal) and non-exonic regions (background noise). Generally speaking, you will need to run the following code to generate observed read count and background read count. 

```
htseq-count [options] <alignment_file> <gtf_file> > Observed_count.txt
htseq-count [options] <alignment_file> <gtf_file_bg> > background_count.txt
```

Details regarding how to use HTSeq can be found here:
http://www-huber.embl.de/HTSeq/doc/count.html

The gtf file used to measure background noise can be downloaded in the gtf folder. 

### XBSeq testing for DE 

After HTSeq procedure, the we will have two measurements for each gene, the observed signal and background noise. The differential expression analysis will be carried out as follows:

```r
# Observe and background are the output matrix from HTSeq
Signal <- estimateRealcount(observe, background)
# conditions are the design matrix for the experiment
XB <- newXBSeqDataSet(Signal,conditions)
XB <- estimateSizeFactorsXBSeq( XB )
XB <-estimateSCV( XB, observe, background, method='pooled', sharingMode='maximum', fitType='local' )
Teststas <- XBSeqTest( XB, levels(conditions)[1L], levels(conditions)[2L], pvals_only=pvals_only )

# Alternatively, all the codes above can be done with a wrapper function XBSeq
Teststats <- XBSeq( observe, background, conditions, method='pooled', sharingMode='maximum', fitType='local', pvals_only=FALSE )
```
# Bug reports
Report bugs as issues on our [GitHub repository](https://github.com/Liuy12/XBSeq) or you can report directly to my email: liuy12@uthscsa.edu.

# Session information 
```r
sessionInfo()
```
```
# R version 3.1.2 (2014-10-31)
# Platform: x86_64-w64-mingw32/x64 (64-bit)

# locale:
# [1] LC_COLLATE=English_United States.1252 
# [2] LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets 
# [6] methods   base     

# other attached packages:
# [1] XBSeq_0.99.0 MASS_7.3-35 

# loaded via a namespace (and not attached):
# [1] Biobase_2.26.0      BiocGenerics_0.12.1
# [3] colorspace_1.2-4    digest_0.6.8       
# [5] ggplot2_1.0.0       grid_3.1.2         
# [7] gtable_0.1.2        lattice_0.20-29    
# [9] locfit_1.5-9.1      matrixStats_0.12.2 
# [11] munsell_0.4.2       parallel_3.1.2     
# [13] plyr_1.8.1          pracma_1.7.9       
# [15] proto_0.3-10        R.methodsS3_1.6.1  
# [17] Rcpp_0.11.3         reshape2_1.4.1     
# [19] scales_0.2.4        stringr_0.6.2      
# [21] tools_3.1.2          
```
# Acknowledgements 
XBSeq is implemented in R based on the source code from DESeq. 
