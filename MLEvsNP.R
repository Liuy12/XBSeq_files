setwd('Dropbox/research/PoissonNB/')
source('Simulation/robust_simulation.R')
Observe<-read.table('Sarkar/allSamples_Genes_ReadCount1.txt',row.names=1,header=T,sep='\t')
Intron<-read.table('Sarkar/allSamples_Genes_shift_ReadCount.txt',row.names=1,header=T,sep='\t')
Observe<-as.matrix(Observe[-c((nrow(Observe)-4):nrow(Observe)),7:9])
grp <- as.factor(rep(0:1,each=10))
data_10v10_1.5fd_0.1de <- NBsim(foldDiff=1.5,dataset=Observe,nTags=5000,group=grp,verbose=T,add.outlier=F,drop.extreme.dispersion=0.1,seed=1990,pDiff=0.1)
index_sample <- data_10v10_1.5fd_0.1de$dataset$dataset.index[data_10v10_1.5fd_0.1de$id_r]
set.seed(1990)
lambda <- rowMeans(Intron[index_sample,])

plusfactor <- 20
Intron_simu_10v10 <- matrix(round(100*rnorm(20*5000,rpois(20*5000,lambda+plusfactor),3)),nrow=5000,ncol=20)
Intron_simu_10v10[Intron_simu_10v10<0] <- 0

Intron_simu_10v10_1 <- matrix(round(rnorm(20*5000, mean=Intron_simu_10v10, sd=0.1*Intron_simu_10v10)),nrow=5000,ncol=20)
Intron_simu_10v10_1[Intron_simu_10v10_1<0]<-0


data_noise_10v10 <- data_10v10_1.5fd_0.1de 
data_noise_10v10$counts <- data_10v10_1.5fd_0.1de$counts + Intron_simu_10v10

pvalue_DESeq <- pval(data_noise_10v10, method='DESeq_pool', count.type='counts', mc.cores=1)
pvalue_XBSeq_noise_10v10_NP <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor( rep( c('C1','C2'), each=10 ) ), fitType='local',  paraMethod = 'NP')
pvalue_XBSeq_noise_10v10_MLE <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor( rep( c('C1','C2'), each=10 ) ), fitType='local',  paraMethod = 'MLE')

pvalue <- pvalue_DESeq
pvalue$pval$XBSeq_NP <- pvalue_XBSeq_noise_10v10_NP$pval
pvalue$pval$XBSeq_MLE <- pvalue_XBSeq_noise_10v10_MLE$pval
pvalue$pval$DESeq <- pvalue_DESeq$pval[[1]]
pvalue$pval$DESeq_pool<-c()

pvalue$method<-c('XBSeq_NP','XBSeq_MLE', 'DESeq')
pvalue$methodVersion<-pvalue$method
pvalue$main<-'10v10_0.1'

roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
fdPlot(pvalue,pout='pval',col=rep(c("black", "blue", "purple"),times=1),lty=rep(1,each=3))
getPower(pvalue, byAveLogCPM = F,pout='pval',col=rep(c("black", "blue", "purple"),times=1),plot=T,lty=rep(1,each=3))

############################# Try baseline 10 replicates
Intron_simu_10v10_baseline <- matrix(rpois(20*5000,lambda), nrow=5000, ncol=20)
Intron_simu_10v10_baseline[Intron_simu_10v10_baseline < 0] <- 0

Intron_simu_10v10_baseline_1 <- matrix(round(rnorm(20*5000, mean=Intron_simu_10v10_baseline, sd=0.1*Intron_simu_10v10_baseline)), nrow=5000, ncol=20)
Intron_simu_10v10_baseline_1[Intron_simu_10v10_baseline_1 < 0] <- 0

data_noise_10v10_baseline <- data_10v10_1.5fd_0.1de 
data_noise_10v10_baseline$counts <- data_10v10_1.5fd_0.1de$counts + Intron_simu_10v10_baseline

pvalue_DESeq_baseline <- pval(data_noise_10v10_baseline, method='DESeq_pool', count.type='counts', mc.cores=1)
pvalue_XBSeq_noise_10v10_NP_baseline <- XBSeq(data_noise_10v10_baseline$counts, Intron_simu_10v10_baseline_1, as.factor( rep( c('C1','C2'), each=10 ) ), fitType='local',  paraMethod = 'NP')
pvalue_XBSeq_noise_10v10_MLE_baseline <- XBSeq(data_noise_10v10_baseline$counts, Intron_simu_10v10_baseline_1, as.factor( rep( c('C1','C2'), each=10 ) ), fitType='local',  paraMethod = 'MLE')

pvalue1 <- pvalue_DESeq_baseline
pvalue1$pval$XBSeq_NP <- pvalue_XBSeq_noise_10v10_NP_baseline$pval
pvalue1$pval$XBSeq_MLE <- pvalue_XBSeq_noise_10v10_MLE_baseline$pval
pvalue1$pval$DESeq <- pvalue_DESeq_baseline$pval[[1]]
pvalue1$pval$DESeq_pool<-c()

pvalue1$method<-c('XBSeq_NP','XBSeq_MLE', 'DESeq')
pvalue1$methodVersion<-pvalue$method
pvalue1$main<-'10v10_0.1'

roPlot(pvalue1,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
fdPlot(pvalue1,pout='pval',col=rep(c("black", "blue", "purple"),times=1),lty=rep(1,each=3))
getPower(pvalue1, byAveLogCPM = F,pout='pval',col=rep(c("black", "blue", "purple"),times=1),plot=T,lty=rep(1,each=3))
