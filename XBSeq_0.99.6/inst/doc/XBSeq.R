## ----eval=FALSE----------------------------------------------------------
#  source('http://www.bioconductor.org/biocLite.R')
#  biocLite("XBSeq")

## ----message = FALSE,eval=FALSE------------------------------------------
#  library("XBSeq")

## ----eval=FALSE----------------------------------------------------------
#  data(ExampleData)

## ----eval=FALSE----------------------------------------------------------
#  head(Observed)
#  head(Background)

## ----tidy=TRUE,eval=FALSE------------------------------------------------
#  conditions <- factor(c(rep('C',3), rep('T', 3)))
#  XB <- XBSeqDataSet(Observed, Background, conditions)

## ---- tidy=TRUE,eval=FALSE-----------------------------------------------
#  XB <- estimateRealCount(XB)
#  XB <- estimateSizeFactors(XB)
#  XB <-estimateSCV( XB, method='pooled', sharingMode='maximum', fitType='local' )

## ----fig.width=5,fig.height=5,eval=FALSE---------------------------------
#  plotSCVEsts(XB)

## ----eval=FALSE----------------------------------------------------------
#  Teststas <- XBSeqTest( XB, levels(conditions)[1L], levels(conditions)[2L] )

## ----fig.width=5,fig.height=5,eval=FALSE---------------------------------
#  MAplot(Teststas, padj = FALSE, pcuff = 0.01, lfccuff = 1)

## ----eval=FALSE,tidy=TRUE------------------------------------------------
#  # Alternatively, all the codes above can be done with a wrapper function XBSeq
#  Teststats <- XBSeq( Observed, Background, conditions, method='pooled', sharingMode='maximum',
#    fitType='local', pvals_only=FALSE )

## ----eval=FALSE----------------------------------------------------------
#  biocLite("DESeq")

## ----message=FALSE,eval=FALSE--------------------------------------------
#  library('DESeq')
#  library('ggplot2')
#  de <- newCountDataSet(Observed, conditions)
#  de <- estimateSizeFactors(de)
#  de <- estimateDispersions(de, method = "pooled", fitType="local")
#  res <- nbinomTest(de, levels(conditions)[1], levels(conditions)[2])

## ----warning=FALSE,message=FALSE,tidy=TRUE, fig.width=5,fig.height=5,eval=FALSE----
#  DE_index_DESeq <- with(res, which(pval<0.01 & abs(log2FoldChange)>1))
#  DE_index_XBSeq <- with(Teststas, which(pval<0.01 & abs(log2FoldChange)>1))
#  DE_index_inters <- intersect(DE_index_DESeq, DE_index_XBSeq)
#  DE_index_DESeq_uniq <- setdiff(DE_index_DESeq, DE_index_XBSeq)
#  DE_plot <- MAplot(Teststas, padj = FALSE, pcuff = 0.01, lfccuff = 1, shape=16)
#  DE_plot + geom_point( data=Teststas[DE_index_inters,], aes(x=baseMean, y=log2FoldChange),
#                        color= 'green', shape=16 ) +
#    geom_point( data=Teststas[DE_index_DESeq_uniq,], aes( x=baseMean, y=log2FoldChange ),
#                color= 'blue', shape=16 )

## ----eval=FALSE----------------------------------------------------------
#  sessionInfo()

