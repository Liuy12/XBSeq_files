# Instruction for how to use the scripts to construct gtf annotation for background noise

Fistly, download refFlat table from UCSC database (http://genome.ucsc.edu) for your organism of interest and create the preliminary list of gene-free regions by using exonFreeRegionShift.pl:

```
 exonFreeRegionShift.pl <-EX exon-GTF file > <-FR gene free region>
```

Then, download tables of (a) all_mrna; (b) ensGene; (c) pseudoYale60Gene; (d) vegaGene;, (e)xenoMrna, and (f) xenoRefGene from UCSC database for your organism of interest and remove regions appear in any of them from the gene-free regions by using GEFRshift.pl:

```
GEFRshift.pl <-G gene-GTF.gtf > <-I intronRegion.tsv> <-T integenicRegion.tsv>
       optional: -m mRNA.bed -x xenoMrna.bed -z xenoRefGene.bed -e ensGene.bed -p pseudoGene.bed -v vegaGene.bed -b 
```

After the two steps, you will be able to generate the background annotation file. 
