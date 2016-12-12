setwd('Dropbox/research/PoissonNB/')
source('Simulation/robust_simulation.R')
Observe<-read.table('Sarkar/allSamples_Genes_ReadCount1.txt',row.names=1,header=T,sep='\t')
Intron<-read.table('Sarkar/allSamples_Genes_shift_ReadCount.txt',row.names=1,header=T,sep='\t')
Observe<-as.matrix(Observe[-c((nrow(Observe)-4):nrow(Observe)),7:9])

##############################################################################
######## examine MLE vs NP
# 2 replicates
pvalue <- Simu_fun(Observe, as.factor(rep(0:1,each=2)), 1.5, 5000 ,1990)
roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
# 0.6133764, 0.6781787, 0.6639529

# 5 replicates
pvalue <- Simu_fun(Observe, as.factor(rep(0:1,each=5)), 1.5, 5000 ,1990)
roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
# 0.6918067, 0.7904809, 0.7855071

# 10 replicates
pvalue <- Simu_fun(Observe, as.factor(rep(0:1,each=10)), 1.5, 5000 ,1990)
roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
# 0.7660507, 0.8393582, 0.839628

# 20 replicates
pvalue <- Simu_fun(Observe, as.factor(rep(0:1,each=20)), 1.5, 5000 ,1990)
roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
# 0.8077591, 0.9038551, 0.9038982

# 50 replicates
pvalue <- Simu_fun(Observe, as.factor(rep(0:1,each=50)), 1.5, 5000 ,1990)
roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
# 0.8568991, 0.9453244, 0.9453244

# 100 replicates
pvalue <- Simu_fun(Observe, as.factor(rep(0:1,each=100)), 1.5, 5000 ,1990)
roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,lty=rep(1,each=3),col=c("black", "blue", "purple"), psi=5)
# 0.9249884, 0.9678244, 0.9678244

################################################################################
##### examine the effect of big_count
fd <- 1.5
ngenes <- 5000
grp <- rep(c(1,2), each = 3)
seed <- 1990

data_10v10_1.5fd_0.1de <- NBsim(foldDiff=fd,dataset=Observe,nTags=ngenes,group=grp,verbose=T,add.outlier=F,drop.extreme.dispersion=0.1,seed=seed,pDiff=0.1, simfunc = 'NB')
index_sample <- data_10v10_1.5fd_0.1de$dataset$dataset.index[data_10v10_1.5fd_0.1de$id_r]
set.seed(seed)
lambda <- rowMeans(Intron[index_sample,])

plusfactor <- 20
ln <- length(grp)
Intron_simu_10v10 <- matrix(round(100*rnorm(ln*ngenes,rpois(ln*ngenes,lambda+plusfactor),3)),nrow=ngenes,ncol=ln)
Intron_simu_10v10[Intron_simu_10v10<0] <- 0

Intron_simu_10v10_1 <- matrix(round(rnorm(ln*ngenes, mean=Intron_simu_10v10, sd=0.1*Intron_simu_10v10)),nrow=ngenes,ncol=ln)
Intron_simu_10v10_1[Intron_simu_10v10_1<0]<-0

data_noise_10v10 <- data_10v10_1.5fd_0.1de 
data_noise_10v10$counts <- data_10v10_1.5fd_0.1de$counts + Intron_simu_10v10

pvalue_XBSeq_noise_10v10 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = max(data_noise_10v10$counts))
pvalue_XBSeq_noise_10v10_bigcount0 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = 0)
pvalue_XBSeq_noise_10v10_bigcount100 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = 100)
pvalue_XBSeq_noise_10v10_bigcount500 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = 500)
pvalue_XBSeq_noise_10v10_bigcount1000 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = 1000)
pvalue_XBSeq_noise_10v10_bigcount5000 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = 5000)
pvalue_XBSeq_noise_10v10_bigcount10000 <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP', big_count = 10000)

pvalue_DESeq2 <- pval(data_noise_10v10, method='DESeq2', count.type='counts', mc.cores=1)
pvalue <- pvalue_DESeq2
pvalue$pval$XBSeq2 <- pvalue_XBSeq_noise_10v10$pval
pvalue$pval$XBSeq2_0 <- pvalue_XBSeq_noise_10v10_bigcount0$pval
pvalue$pval$XBSeq2_100 <- pvalue_XBSeq_noise_10v10_bigcount100$pval
pvalue$pval$XBSeq2_500 <- pvalue_XBSeq_noise_10v10_bigcount500$pval
pvalue$pval$XBSeq2_1000 <- pvalue_XBSeq_noise_10v10_bigcount1000$pval
pvalue$pval$XBSeq2_5000 <- pvalue_XBSeq_noise_10v10_bigcount5000$pval
pvalue$pval$XBSeq2_10000 <- pvalue_XBSeq_noise_10v10_bigcount10000$pval


roPlot(pvalue,pout='pval',threshold=0.05,plot.max.fpr=1,col=1:8, psi=5)

pvalue$method<-1:8
pvalue$methodVersion<-pvalue$method
pvalue$main<-'10v10_0.1'




Simu_fun <- function(data, grp, fd, ngenes, seed){
  data_10v10_1.5fd_0.1de <- NBsim(foldDiff=fd,dataset=data,nTags=ngenes,group=grp,verbose=T,add.outlier=F,drop.extreme.dispersion=0.1,seed=seed,pDiff=0.1, simfunc = 'NB')
  index_sample <- data_10v10_1.5fd_0.1de$dataset$dataset.index[data_10v10_1.5fd_0.1de$id_r]
  set.seed(seed)
  lambda <- rowMeans(Intron[index_sample,])
  
  plusfactor <- 20
  ln <- length(grp)
  Intron_simu_10v10 <- matrix(round(100*rnorm(ln*ngenes,rpois(ln*ngenes,lambda+plusfactor),3)),nrow=ngenes,ncol=ln)
  Intron_simu_10v10[Intron_simu_10v10<0] <- 0
  
  Intron_simu_10v10_1 <- matrix(round(rnorm(ln*ngenes, mean=Intron_simu_10v10, sd=0.1*Intron_simu_10v10)),nrow=ngenes,ncol=ln)
  Intron_simu_10v10_1[Intron_simu_10v10_1<0]<-0
  
  data_noise_10v10 <- data_10v10_1.5fd_0.1de 
  data_noise_10v10$counts <- data_10v10_1.5fd_0.1de$counts + Intron_simu_10v10
  
  cat('Statistical testing using DESeq2', '\n')
  pvalue_DESeq2 <- pval(data_noise_10v10, method='DESeq2', count.type='counts', mc.cores=1)
  cat('Statistical testing using XBSeq NP', '\n')
  pvalue_XBSeq_noise_10v10_NP <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'NP')
  cat('Statistical testing using XBSeq MLE', '\n')
  pvalue_XBSeq_noise_10v10_MLE <- XBSeq(data_noise_10v10$counts, Intron_simu_10v10_1, as.factor(grp), fitType='local',  paraMethod = 'MLE')
  
  pvalue <- pvalue_DESeq2
  pvalue$pval$XBSeq_NP <- pvalue_XBSeq_noise_10v10_NP$pval
  pvalue$pval$XBSeq_MLE <- pvalue_XBSeq_noise_10v10_MLE$pval
  
  pvalue$method<-c('DESeq2', 'XBSeq_NP','XBSeq_MLE')
  pvalue$methodVersion<-pvalue$method
  pvalue$main<-'10v10_0.1'
  return(pvalue)
}




#fdPlot(pvalue,pout='pval',col=rep(c("black", "blue", "purple"),times=1),lty=rep(1,each=3))
#getPower(pvalue, byAveLogCPM = F,pout='pval',col=rep(c("black", "blue", "purple"),times=1),plot=T,lty=rep(1,each=3))
