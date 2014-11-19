# Introduction 

XBSeq is a novel algorithm for testing RNA-seq differential expression (DE), where a statistical model was established
based on the assumption that observed signals are the convolution of true expression signals and sequencing noises. The
mapped reads in non-exonic regions are considered as sequencing noises, which follows a Poisson distribution. Given
measureable observed and noise signals from RNA-seq data, true expression signals, assuming governed by the negative
binomial distribution, can be delineated and thus the accurate detection of differential expressed genes
# Installation 

Currently we are still finalizing the documentation for XBSeq to be released in bioconductor. Right now the algorithm is implemented in R. To use the algorithm, first you need to download the XBSeq.R file into your local directory and then load into your workspace by:

```r
source('your_local_directory/XBSeq')
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

After HTSeq procedure, the 
