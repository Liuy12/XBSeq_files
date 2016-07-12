# gtf Files for XBSeq
There are two gtf files for each genome build

* gtf file ends with gene_shift.gtf
This file is for use to generate background noise for a RNA-seq experiment by using HTSeq

* gtf file ends with .gtf
This file is for use to generate observed signal for a RNA-seq experiment by using HTSeq

In order to use XBSeq for testing DE, after sequence alignment, we need to run HTSeq twice to measure the reads mapped to exonic regions (observed signal) and non-exonic regions (background noise). Generally speaking, you will need to run the following code to generate observed signal and background noise. 

```{r,engine='python',eval=FALSE}
htseq-count [options] <alignment_file> <gtf_file> > Observed_count.txt
htseq-count [options] <alignment_file> <gtf_file_bg> > background_count.txt
```
Details regarding how HTSeq works can be found here: http://www-huber.embl.de/HTSeq/doc/count.html

For more details please refer to XBSeq package: https://bioconductor.org/packages/devel/bioc/html/XBSeq.html

