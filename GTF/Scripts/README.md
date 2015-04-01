# Instruction for how to use the scripts to construct gtf annotation for background noise

Fistly, download refFlat table from UCSC database (http://genome.ucsc.edu) and create the preliminary list of gene-free regions 
by using exonFreeRegionShift.pl

```r
 exonFreeRegionShift.pl <-EX exon-GTF file > <-FR gene free region>
```
	 \item {Download tables of (a) all_mrna; (b) ensGene; (c) pseudoYale60Gene; (d) vegaGene;, (e)xenoMrna, and (f) xenoRefGene from UCSC database and remove 
	 regions appear in any of them from the gene-free regions,
