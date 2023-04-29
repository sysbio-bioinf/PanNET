
library("BoolNet") #we used version 2.1.7
library("qgraph") #we used version 1.9.3
library("igraph") #we used version 1.3.5
library("pROC") #we used version 1.18.0
library("ggplot2") #we used version 3.4.1
library("ggpubr") #we used version 0.6.0
library("poweRlaw") #we used version 0.70.6


####Figure 6 
GSE73514NormalTumor<-read.csv("GSE73514StagesNormalTumor.csv", sep = ",", header = T)
as.data.frame(GSE73514NormalTumor)
genesonemore<-colnames(GSE73514NormalTumor)
genes<-genesonemore[-1]

#Routine for plotting the box plots
pdf("Results.pdf")
pvalvec <- rep(NA, ncol(GSE73514NormalTumor)-1)
for (a in genes){
  #wilcoxon testing
  print("Column nr being used:")
  print(which(genes==a)+1)
  pvalvec[which(genes==a)]<-wilcox.test(GSE73514NormalTumor[,which(genes==a)+1]~GSE73514NormalTumor$Stage,alternative=c("two.sided"), paired=F, conf.level=0.95)$p.value
  print(test)
  
  #Boxplots
  g<-ggplot(data= subset(GSE73514NormalTumor , select = c("Stage",a)), aes_string(x="Stage", y=a)) + geom_boxplot(fill="#A6C851") + geom_jitter() + ggtitle(paste(a, "p-value=",pvalvec[which(genes==a)]))#+ labs(x="",y=paste(a, "expression") + geom_jitter() #+ geom_hline(yintercept = 8.4, linetype="dashed", color="#AB1212")
  print(g)
  
  #Roc curves
  roc<-  roc(GSE73514NormalTumor$Stage, GSE73514NormalTumor[,which(genes==a)+1], percent=TRUE,partial.auc.correct=FALSE,
      plot=TRUE,thresholds="best", 
      print.thres="best", print.auc = TRUE, main=a)
  print(roc)
}
dev.off()

pvalueSummary <- cbind(genes, pvalvec)

