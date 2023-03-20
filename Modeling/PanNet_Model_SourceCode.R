####### Source code ####### 
#### Werle et al. 2023 - A systems biology approach to define mechanisms, phenotypes, and drivers in PanNETs with a personalized perspective #####

#Set working directory
setwd("~/User/SourceCode/")

#Load required packages
#We run the simulations with R version 4.2.2 (2022-10-31)
library("BoolNet") #we used version 2.1.7
library("qgraph") #we used version 1.9.3
library("igraph") #we used version 1.3.5
library("pROC") #we used version 1.18.0
library("ggplot2") #we used version 3.4.1
library("ggpubr") #we used version 0.6.0
library("poweRlaw") #we used version 0.70.6

#Load network
PanNET<-loadNetwork("PanNet_Model.txt", symbolic = T) ###network with time delays
PanNET_noDelays<-loadNetwork("PanNet_Model_NoDelays.txt") ###network without time delays
####To reduce the number of nodes in the network, linear connections were removed by time-delays. 
####The network with no time delays included these nodes and thus has a slightly higher number of nodes
#56 genes in reduced version, 66 in NoDelays version

#Export to SBML
toSBML(PanNET_noDelays, file="PanNETnoDelay.sbml")

#Adjust size for attractor plotting
par(xpd = T)
par(mar = c(1,7,7,1))

###################################################################################################################
####Main
###################################################################################################################
####Figure 3b
######The code for the basin estimation on which these pie charts are based can be found in the Supplement part Figure S4-S7

#WT
slices<- c(99, 1)
lbls<- c("Quiescence", "Angiogenesis, Detachment, Proliferation") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#8A2BE2"), main = "WT PanNET")

#DAXX loss
slices<- c(98.8, 1.1, 0.9)
lbls<- c("Quiescence", "G0-alert, detachment", "Angiogenesis, Detachment, Proliferation") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#D6D6D6", "#8A2BE2"), main = "DAXX loss")

#TSC loss
slices<- c(99, 1)
lbls<- c("G0-Alert", "Angiogenesis, Detachment, Proliferation")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#2E7DFF", "#8A2BE2"), main = "TSC loss")

#MEN1 loss
slices<- c(21.68, 19.93, 13.05, 45.33)
lbls<- c("Quiescence", "G0-alert","Angionesis, GO-alert", "Angiogenesis, Detachment, Proliferation") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#2E7DFF", "#DA8C85", "#8A2BE2"), main = "MEN-1 loss")

###################################################################################################################
####Figure 4c (in addition to the pie charts already in Figure 3b)
######The code for the basin estimation on which this pie charts is based can be found in the Supplement part Figure S9

#MEN1 and DAXX loss
slices<- c(21.46, 20.14, 46.12, 12.26 )
lbls<- c("Quiescence", "G0-alert", "Angiogenesis, Detachment, Proliferation", "Angiogenesis, G0-alert") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#2E7DFF", "#8A2BE2", "#DA8C85"), main = "DAXX and MEN1 loss")

###################################################################################################################
####Figure 5b
slices<- c(14.3, 67.9, 17.9)
lbls<- c("Constant activity", "Activity change", "Partial activity change")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("cyan4", "hotpink4", "orange"), main = "Comparison node activity WT quiescence vs severe")

####Figure 5c
slices<- c(69.64, 30.36)
lbls<- c("validated", "not validated")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("lightgoldenrod2", "paleturquoise3"), main = "Validation nodes")

###################################################################################################################
####Figure 6 
GSE73514NormalTumor<-read.csv("GSE73514StagesNormalTumor.csv", sep = ",", header = T)

#CCNE1
rocCCNE1 <- roc(GSE73514NormalTumor$Stage, GSE73514NormalTumor$Ccne1, percent=TRUE,partial.auc.correct=FALSE,
                plot=TRUE,thresholds="best", 
                print.thres="best", print.auc = TRUE, main="Ccne1")
wilcox.test(GSE73514NormalTumor$Ccne1~GSE73514NormalTumor$Stage,alterantive=c("two-sided"), paired=F, conf.level=0.95)

Plot<-ggplot(GSE73514NormalTumor, aes(Stage, Ccne1))
Plot + geom_boxplot(fill="#A6C851")  + labs(x="",y="Ccne1 expression")   +geom_jitter() +geom_hline(yintercept = 8.4, linetype="dashed", color="#AB1212")

#E2F
rocE2F <- roc(GSE73514NormalTumor$Stage, GSE73514NormalTumor$E2f1, percent=TRUE,partial.auc.correct=FALSE,
              plot=TRUE,thresholds="best", 
              print.thres="best", print.auc = TRUE, main="E2f1")

wilcox.test(GSE73514NormalTumor$E2f1~GSE73514NormalTumor$Stage,alterantive=c("two-sided"), paired=F, conf.level=0.95)

Plot<-ggplot(GSE73514NormalTumor, aes(Stage, E2f1))
Plot + geom_boxplot(fill="#A6C851")  + labs(x="",y="E2f1 expression")   +geom_jitter() +geom_hline(yintercept = 8.7, linetype="dashed", color="#AB1212")

#CDH1
rocCDH1 <- roc(GSE73514NormalTumor$Stage, GSE73514NormalTumor$Cdh1, percent=TRUE,partial.auc.correct=FALSE,
               plot=TRUE,thresholds="best", 
               print.thres="all", print.auc = TRUE, main="Cdh1")
wilcox.test(GSE73514NormalTumor$Cdh1~GSE73514NormalTumor$Stage,alterantive=c("two-sided"), paired=F, conf.level=0.95)

Plot<-ggplot(GSE73514NormalTumor, aes(Stage, Cdh1))
Plot + geom_boxplot(fill="#A6C851")  + labs(x="",y="Cdh1 expression")   +geom_jitter() +geom_hline(yintercept = 10.4, linetype="dashed", color="#AB1212")

#CIP2A
rocCIP2A <- roc(GSE73514NormalTumor$Stage, GSE73514NormalTumor$Cip2a, percent=TRUE,partial.auc.correct=FALSE,
                plot=TRUE,thresholds="best", 
                print.thres="best", print.auc = TRUE, main="Cip2a")

wilcox.test(GSE73514NormalTumor$Cip2a~GSE73514NormalTumor$Stage,alterantive=c("two-sided"), paired=F, conf.level=0.95)

Plot<-ggplot(GSE73514NormalTumor, aes(Stage, Cip2a))
Plot + geom_boxplot(fill="#A6C851")  + labs(x="",y="Cip2a expression")   +geom_jitter() +geom_hline(yintercept = 8.1, linetype="dashed", color="#AB1212")

#IGF
rocIGF <- roc(GSE73514NormalTumor$Stage, GSE73514NormalTumor$Igf1, percent=TRUE,partial.auc.correct=FALSE,
              plot=TRUE,thresholds="best", 
              print.thres="best", print.auc = TRUE, main="Igf1")

wilcox.test(GSE73514NormalTumor$Igf1~GSE73514NormalTumor$Stage,alterantive=c("two-sided"), paired=F, conf.level=0.95)

Plot<-ggplot(GSE73514NormalTumor, aes(Stage, Igf1))
Plot + geom_boxplot(fill="#A6C851")  + labs(x="",y="Igf1 expression")   +geom_jitter() +geom_hline(yintercept = 7.3, linetype="dashed", color="#AB1212")

###################################################################################################################
#Figure 7

#########This part needs to be run on a server. We provide the output files for 100 million runs as RDS files below. Thus, this part can be skipped
N <- 100000000
netfile <- "PanNet_Model_NoDelays.txt"
SBML <- loadNetwork(netfile)
CCNE1index <- which(SBML$genes == "CCNE1")
E2Findex <- which(SBML$genes == "E2F")
TJindex <- which(SBML$genes == "TJ")
driverGenes <- rep(0, length(SBML$genes))
cancerStateCount <- 0
names(driverGenes) <- SBML$genes
MENlossSBML <- SBML
MENlossSBML$fixed[49] <- 0

#Menloss <- TRUE #to simulate MEN1 loss
Menloss <- FALSE #to simulate the WT

for (s in 1:N){
  if (s %% 100 == 0){
    print(paste0("Run: ", s, "/", N))
    print(paste0("The ratio of cancer states so far is ", cancerStateCount, "/", s, " = ", round(cancerStateCount/s,3)))
  }
  if (Menloss){net <- MENlossSBML} else {net <- SBML}
  startState <- sample(c(0,1), size=length(net$genes), replace = TRUE)
  if (Menloss){startState[49] <- 0}
  attr <- getAttractors(net, method="chosen", startStates = list(startState))
  binaryAttr <- t(as.matrix(getAttractorSequence(attr, 1)))
  attrLength <- ncol(binaryAttr)
  if (sum(binaryAttr[CCNE1index,]) == attrLength & 
      sum(binaryAttr[E2Findex,]) == attrLength &
      sum(binaryAttr[TJindex,]) < attrLength){
    cancerStateCount <- cancerStateCount + 1
    driverGenes <- driverGenes + startState
  }
}

#select the RDS file you prefer to save depending on the simulation setup --> uncomment the following rows to do it
#saveRDS(driverGenes/cancerStateCount, file="MEN1KOprobs_100000000states_new.RDS")
#saveRDS(driverGenes/cancerStateCount, file="WTprobs_100000000states_new.RDS")

#########End part that can be skipped

#To retrieve the graphics of the manuscript the code can be run from here on just loading the RDS files
WT_probs <- readRDS("WTprobs_100000000states.RDS")
MEN1KO_probs <- readRDS("MEN1KOprobs_100000000states.RDS")

networkNameOrder <- names(WT_probs)
sd_WT <- sd(WT_probs)
sd_WT
sd_WT_fromBaseline <- sd(abs(WT_probs-0.5))
sd_WT_fromBaseline
sd_MEN1KO <- sd(MEN1KO_probs)
sd_MEN1KO
sd_MEN1KO_fromBaseline <- sd(abs(MEN1KO_probs-0.5))
sd_MEN1KO_fromBaseline

barplot(WT_probs, ylim = c(0,1))
abline(a=0.5, b=0, col="red")

#WT
df <- data.frame(
  name=names(WT_probs),
  value=WT_probs
)

df <- data.frame(
  name=names(WT_probs),
  value=WT_probs#driverGenes/cancerStateCount
)
df$name <- factor(df$name, levels = df$name)

ggplot(df, aes(x=name, y=value)) + ylim(c(0,1)) + geom_hline(yintercept=0.5, col="red", size=1) + 
  xlab("Gene") + ylab("Probability of gene being 1 in start state when entering cancer attractor") +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept=0.5+sd_WT_fromBaseline, col="blue") + geom_hline(yintercept=0.5-sd_WT_fromBaseline, col="blue") #+ 


#MEN-1 loss
df <- data.frame(
  name=names(MEN1KO_probs),
  value=MEN1KO_probs
)

df <- data.frame(
  name=names(MEN1KO_probs),
  value=MEN1KO_probs#driverGenes/cancerStateCount
)
df$name <- factor(df$name, levels = df$name)

ggplot(df, aes(x=name, y=value)) + ylim(c(0,1)) + geom_hline(yintercept=0.5, col="red", size=1) + 
  xlab("Gene") + ylab("Probability of gene being 1 in start state when entering cancer attractor") +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept=0.5+sd_MEN1KO_fromBaseline, col="blue") + geom_hline(yintercept=0.5-sd_MEN1KO_fromBaseline, col="blue") #+ 


###################################################################################################################
#Figure 8a
WT<-read.csv("PhenotypicalLandscapeWT.csv", sep = ",", header = T)
WT$Phenotype<-factor(WT$Phenotype, levels = c("Quiescence", "GO-alert", "Proliferation", "Detachement", "Angiogenesis"))
WT$Driver<-factor(WT$Driver, levels = c("WT", "AKT", "CIP2A", "ERK", "FOXO1", "FOXO3", "GLI1", "GSK3B", "MEN1", "mTORC2", "RSK", "PP2A", "TP53"))

ggplot(data=WT, aes(x=Driver, y=Distribution, fill=Phenotype)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("#96C217", "#2E7DFF", "#CFA7FF","#F5D586", "#DA8C85")) +
  geom_line(aes(x = Driver, y = 10*NumberAttractors), size = 1, color="#B467BA", group = 1) + 
  geom_line(aes(x = Driver, y = 10*NumberPhenotpes), size = 1, color="#989898", group = 1) + 
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Number of attractors and phenotypes")) 

#Figure 8b
MEN1<-read.csv("Men1AllIn.csv", sep = ",", header = T)
MEN1$Phenotype<-factor(MEN1$Phenotype, levels = c("Quiescence", "GO-alert", "Proliferation", "Detachement", "Angiogenesis"))
MEN1$Driver<-factor(MEN1$Driver, levels = c("MEN1", "AKT", "CIP2A", "E2F", "ETS1","PP2A", "GLI1", "TP53", "ERK"))

ggplot(data=MEN1, aes(x=Driver, y=Distribution, fill=Phenotype)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("#96C217", "#2E7DFF", "#CFA7FF","#F5D586", "#DA8C85")) +
  geom_line(aes(x = Driver, y = 10*Number_attractors), size = 1, color="#B467BA", group = 1) + 
  geom_line(aes(x = Driver, y = 10*Number_phenotpes), size = 1, color="#989898", group = 1) + 
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Number of attractors and phenotypes")) 

#######Analysis############
#Figure 8a
AKToe<-fixGenes(PanNET, "AKT", 1)
attrWTAKToe<-simulateSymbolicModel(AKToe, startStates = 10000)
plotAttractors(attrWTAKToe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with AKT overexpression", drawLegend = F , allInOnePlot = T)

CIP2Aoe<-fixGenes(PanNET, "CIP2A", 1)
attrWTCIP2Aoe<-simulateSymbolicModel(CIP2Aoe, startStates = 10000)
plotAttractors(attrWTCIP2Aoe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with CIP2A overexpression", drawLegend = F , allInOnePlot = T)

ERKoe<-fixGenes(PanNET, "ERK", 1)
attrWTERKoe<-simulateSymbolicModel(ERKoe, startStates = 10000)
plotAttractors(attrWTERKoe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with ERK overexpression", drawLegend = F , allInOnePlot = T)

FOXO1ko<-fixGenes(PanNET, "FOXO1", 0)
attrWTFOXO1ko<-simulateSymbolicModel(FOXO1ko, startStates = 10000)
plotAttractors(attrWTFOXO1ko, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with FOXO1 knockout", drawLegend = F , allInOnePlot = T)

FOXO3ko<-fixGenes(PanNET, "FOXO3", 0)
attrWTFOXO3ko<-simulateSymbolicModel(FOXO3ko, startStates = 10000)
plotAttractors(attrWTFOXO3ko, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with FOXO3 knockout", drawLegend = F , allInOnePlot = T)

Gli1oe<-fixGenes(PanNET, "GLI1", 1)
attrWTGli1oe<-simulateSymbolicModel(Gli1oe, startStates = 10000)
plotAttractors(attrWTGli1oe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with GLI1 overexpression", drawLegend = F , allInOnePlot = T)

GSK3Bko<-fixGenes(PanNET, "GSK3B", 0)
attrWTGSK3Bko<-simulateSymbolicModel(GSK3Bko, startStates = 10000)
plotAttractors(attrWTGSK3Bko, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with GSK3B knockout", drawLegend = F , allInOnePlot = T)

mTORC2oe<-fixGenes(PanNET, "mTORC2", 1)
attrWTmTORC2oe<-simulateSymbolicModel(mTORC2oe, startStates = 10000)
plotAttractors(attrWTmTORC2oe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with mTORC2 overexpression", drawLegend = F , allInOnePlot = T)

PP2Ako<-fixGenes(PanNET, "PP2A", 0)
attrWTPP2Ako<-simulateSymbolicModel(PP2Ako, startStates = 10000)
plotAttractors(attrWTPP2Ako, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with PP2A knockout", drawLegend = F , allInOnePlot = T)

TP53ko<-fixGenes(PanNET, "TP53", 0)
attrWTTP53ko<-simulateSymbolicModel(TP53ko, startStates = 10000)
plotAttractors(attrWTTP53ko, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with TP53 knockout", drawLegend = F , allInOnePlot = T)

RSKko<-fixGenes(PanNET, "RSK", 0)
attrWTRSKko<-simulateSymbolicModel(RSKko, startStates = 10000)
plotAttractors(attrWTRSKko, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors WT with RSK Knockout", drawLegend = F , allInOnePlot = T)


#Figure 8b
MEN1loss<-fixGenes(PanNET, "MEN1", 0)
MEN1lossAKToe<-fixGenes(MEN1loss, "AKT", 1)
attrMEN1lossAKToe<-getAttractors(MEN1lossAKToe, startStates = 10000)
plotAttractors(attrMEN1lossAKToe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss AKT overexpression", drawLegend = F , allInOnePlot = T)

MEN1lossPP2Ako<-fixGenes(MEN1loss, "PP2A", 0)
attrMEN1lossPP2Ako<-getAttractors(MEN1lossPP2Ako, startStates = 10000)
plotAttractors(attrMEN1lossPP2Ako, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss PP2A knockdown", drawLegend = F , allInOnePlot = T)

MEN1lossCIP2Aoe<-fixGenes(MEN1loss, "CIP2A", 1)
attrMEN1lossCIP2Aoe<-getAttractors(MEN1lossCIP2Aoe, startStates = 10000)
plotAttractors(attrMEN1lossCIP2Aoe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss CIP2A overexpression", drawLegend = F , allInOnePlot = T)

MEN1lossERKoe<-fixGenes(MEN1loss, "ERK", 1)
attrMEN1lossERKoe<-getAttractors(MEN1lossERKoe, startStates = 10000)
plotAttractors(attrMEN1lossERKoe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss ERK overexpression", drawLegend = F , allInOnePlot = T)

MEN1lossE2Foe<-fixGenes(MEN1loss, "E2F", 1)
attrMEN1lossE2Foe<-getAttractors(MEN1lossE2Foe, startStates = 10000)
plotAttractors(attrMEN1lossE2Foe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss E2F overexpression", drawLegend = F , allInOnePlot = T)

MEN1lossGLI1oe<-fixGenes(MEN1loss, "GLI1", 1)
attrMEN1lossGLI1oe<-getAttractors(MEN1lossGLI1oe, startStates = 10000)
plotAttractors(attrMEN1lossGLI1oe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss GLI1 overexpression", drawLegend = F , allInOnePlot = T)

MEN1lossNoDelay<-fixGenes(PanNET_noDelays, "MEN1", 0)
MEN1lossETS1oe<-fixGenes(MEN1lossNoDelay, "ETS1", 1)
attrMEN1lossETS1oe<-getAttractors(MEN1lossETS1oe, startStates = 10000)
plotAttractors(attrMEN1lossETS1oe, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss ETS1 overexpression", drawLegend = F , allInOnePlot = T)

MEN1lossTP53ko<-fixGenes(MEN1loss, "TP53", 0)
attrMEN1lossTP53ko<-getAttractors(MEN1lossTP53ko, startStates = 10000)
plotAttractors(attrMEN1lossTP53ko, onColor = "olivedrab3", offColor = "slategrey", title = "Attractors MEN1 loss TP53 knockdown", drawLegend = F , allInOnePlot = T)

###################################################################################################################
#Figure 9a
#unperturbed see code pie charts Figure 3b

slices<- c(96.4, 3.6)
lbls<- c("Quiescence", "Detachment, proliferation") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#CFA7FF"), main = "WT mTORC1 KO")

slices<- c(97, 3)
lbls<- c("Quiescence", "Detachment, proliferation") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#CFA7FF"), main = "DAXX mTORC1 KO")

slices<- c(96, 4)
lbls<- c("Quiescence", "Detachment, proliferation")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#CFA7FF"), main = "TSC mTORC1 KO")

slices<- c(69.32, 13.58, 17.1)
lbls<- c("Quiescence", "Detachment", "Detachment, proliferation")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#F5D586", "#CFA7FF"), main = "MEN-1 loss mTORC1 KO")

slices<- c(99.46, 0.54)
lbls<- c("Quiescence", "Detachment, Proliferation")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#CFA7FF"), main = "mTORC1 and mTORC2 KO WT")

slices<- c(99.46, 0.54)
lbls<- c("Quiescence", "Detachment, Proliferation") 
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#CFA7FF"), main = "mTORC1 and mTORC2 KO DAXX")

slices<- c(99.54, 0.46)
lbls<- c("Quiescence", "Detachment, Proliferation")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#CFA7FF"), main = "mTORC1 and mTORC2 KO TSC")


slices<- c(82.28, 16.47, 1.25)
lbls<- c("Quiescence", "Detachment", "Detachment, proliferation")
num <- round(slices/sum(slices)*100)
lbls <- paste(lbls, num)
lbls <- paste(lbls, "%", sep = "")
pie(slices, labels = lbls, col=c("#9ACD32", "#F5D586", "#CFA7FF"), main = "mTORC1 and mTORC2 KO MEN-1")

###################################################################################################################
#Figure 9b

mTORC1ko<-fixGenes(PanNET, "mTORC1", 0)

#WT mTORC1 inhibition to resistance
#WT single-state attractor
WT1prolif1<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)
mTORC1resistance <- plotSequence (mTORC1ko, startState = WT1prolif1, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)

#Sequence cyclic attractor 1 WT
WT1prolif2<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)
WT1prolif3<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SETD7 = 0), default = 1)
WT1prolif4<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolif5<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolif6<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK= 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolif7<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK= 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, HIF1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)
WT1prolif8<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, STAT3 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)

#Cascades from a state of the cyclic WT attractor 1 to resistance
mTORC1resistance2 <- plotSequence (mTORC1ko, startState = WT1prolif2, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance3 <- plotSequence (mTORC1ko, startState = WT1prolif3, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance4 <- plotSequence (mTORC1ko, startState = WT1prolif4, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance5 <- plotSequence (mTORC1ko, startState = WT1prolif5, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance6 <- plotSequence (mTORC1ko, startState = WT1prolif6, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance7 <- plotSequence (mTORC1ko, startState = WT1prolif7, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance8 <- plotSequence (mTORC1ko, startState = WT1prolif8, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
#--> all initial states lead to the same phenotype

#Sequence cyclic attractor 2 WT
WT1prolif2_1<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, STAT3 = 0,  HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)
WT1prolif2_2<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0,ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SETD7 = 0), default = 1)
WT1prolif2_3<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)
WT1prolif2_4<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0,  MDM2 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, CDH1=0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)
WT1prolif2_5<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolif2_6<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, STAT3 = 0,  HIF1 = 0, FOXO1 = 0, FOXO3 = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1= 0), default = 1)
WT1prolif2_7<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0,  MDM2 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SETD7= 0), default = 1)
WT1prolif2_8<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)

#Cascades from a state of the cyclic WT attractor 2 to resistance
mTORC1resistance2_1 <- plotSequence (mTORC1ko, startState = WT1prolif2_1, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_2 <- plotSequence (mTORC1ko, startState = WT1prolif2_2, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_3 <- plotSequence (mTORC1ko, startState = WT1prolif2_3, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_4 <- plotSequence (mTORC1ko, startState = WT1prolif2_4, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_5 <- plotSequence (mTORC1ko, startState = WT1prolif2_5, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_6 <- plotSequence (mTORC1ko, startState = WT1prolif2_6, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_7 <- plotSequence (mTORC1ko, startState = WT1prolif2_7, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistance2_8 <- plotSequence (mTORC1ko, startState = WT1prolif2_8, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
#################################

#Cascade WT mTORC1 and mTORC2 inhibition
mTORC1mTORC2ko<-fixGenes(mTORC1ko, "mTORC2", 0)
WT1prolif1<-generateState(PanNET, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)
mTORC1ANDmTORC2resistance <- plotSequence (mTORC1mTORC2ko, startState = WT1prolif1, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade multi mTORC resistance in WT", drawLegend = F)

#################################
##MEN1 loss mTORC1 inhibition to resistance
mTORC1MEN1ko<-fixGenes(mTORC1ko, "MEN1", 0)
MEN1lossResistant<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0,  SETD7 = 0, PDX1=0), default = 1)
MEN1ANDmTORC1resistance <- plotSequence (mTORC1MEN1ko, startState = MEN1lossResistant, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss PanNETs", drawLegend = F)

#Sequence cyclic attractor 1 MEN1 loss
WT1prolifM2_1<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, STAT3 = 0,  HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)
WT1prolifM2_2<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0,ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SETD7 = 0), default = 1)
WT1prolifM2_3<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)
WT1prolifM2_4<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0,  MDM2 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, CDH1=0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)
WT1prolifM2_5<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolifM2_6<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, STAT3 = 0,  HIF1 = 0, FOXO1 = 0, FOXO3 = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1= 0), default = 1)
WT1prolifM2_7<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0,  MDM2 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SETD7= 0), default = 1)
WT1prolifM2_8<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0,TJ = 0,  RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)

#Cascades from a state of the cyclic attractor 1 MEN1 loss to resistance (same sequence as the 8 states attractor of the WT)
mTORC1resistanceMen_1 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_1, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMen_2 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_2, onColor = "olivedrab3", offColor = "slategrey", title = "CCascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMEN_3 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_3, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMEN_4 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_4, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMEN_5 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_5, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMEN_6 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_6, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMEN_7 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_7, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)
mTORC1resistanceMEN_8 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2_8, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade treatment with mTORC1 inhibitor in MEN1 loss pNETs", drawLegend = F)

#Sequence cyclic attractor 2 MEN1 loss (same sequence as the 8 states attractor of the WT)
WT1prolifM2<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)
WT1prolifM3<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SETD7 = 0), default = 1)
WT1prolifM4<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolifM5<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolifM6<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK= 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0), default = 1)
WT1prolifM7<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK= 0, RSK = 0, TSC = 0, AMPK = 0, TP53 = 0, FourEBP1 = 0, HIF1 = 0, STAT3 = 0, FOXO1 = 0, FOXO3 = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0), default = 1)
WT1prolifM8<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, RSK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, STAT3 = 0, HIF1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, PDX1 = 0, SETD7 = 0), default = 1)

#Cascades from a state of the cyclic attractor 2 MEN1 loss to resistance
mTORC1resistanceM2 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM2, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistanceM3 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM3, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistanceM4 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM4, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistanceM5 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM5, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistanceM6 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM6, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistanceM7 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM7, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)
mTORC1resistanceM8 <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM8, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)

#6-states attractors and 30 states attractors are all oscillations of the previous 8 states attractors and the one reported below 
WT1prolifM1_6states<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PI3K=0, PTEN = 0, AKT=0,  mTORC2 = 0, GSK3B = 0, PP2A = 0, RAS=0,  RAF = 0,  RSK = 0,  AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, HIF1 = 0, SLUG1 = 0,  MYC = 0, TJ= 0, MEN1 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0, SKP2 = 0, FAK1 = 0, PDX1 = 0, SETD7 = 0), default = 1)
mTORC1resistanceM1_6states <- plotSequence (mTORC1MEN1ko, startState = WT1prolifM1_6states, onColor = "olivedrab3", offColor = "slategrey", title = "Cascade mTORC1 resistance in WT", drawLegend = F)

#################################
##MEN1 loss mTORC1 & mTORC2 inhibition to resistance
mTORC1MEN1ko<-fixGenes(mTORC1ko, "MEN1", 0)
mTORC1mTORC2MEN1ko<-fixGenes(mTORC1MEN1ko, "mTORC2", 0)
MEN1lossResistant<-generateState(mTORC1MEN1ko, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0,  SETD7 = 0), default = 1)
MEN1ANDmTORC1resistance <- plotSequence (mTORC1mTORC2MEN1ko, startState = MEN1lossResistant, onColor = "olivedrab3", offColor = "slategrey", title = "pan mTORC inhibitors in MEN1 loss PanNETs", drawLegend = F)

#################################
##MEN1 loss CIP2A inhibition to resistance
MEN1ko<-fixGenes(PanNET, "MEN1", 0)
MEN1CIP2Ako<-fixGenes(MEN1ko, "CIP2A", 0)
MEN1lossResistant<-generateState(MEN1CIP2Ako, specs = c(IRS = 0, PTEN = 0, mTORC2 = 0, GSK3B = 0, PP2A = 0, RAF = 0, MEK = 0, ERK = 0, TSC = 0, AMPK = 0, REDD = 0, TP53 = 0, MDM2 = 0, FourEBP1 = 0, FOXO1 = 0, FOXO3 = 0, CDH1 = 0, TJ = 0, RASSF1A = 0, MEN1 = 0, P27 = 0, P18 = 0, P21 = 0, RB = 0, DAXX = 0,  SETD7 = 0), default = 1)
MEN1withCIP2A <- plotSequence (MEN1CIP2Ako, startState = MEN1lossResistant, onColor = "olivedrab3", offColor = "slategrey", title = "CIP2A inhibitor in MEN1 loss PanNETs", drawLegend = F)

###################################################################################################################
###################################################################################################################
#####Supplement 
###################################################################################################################
###################################################################################################################
##Test for scale-freeness

if(!require("poweRlaw",character.only = TRUE)) install.packages("poweRlaw")
library("poweRlaw",character.only = TRUE)

readNetwork <- function(nwFileName, nwName) {
  inFile <- paste(nwFileName,sep="")
  return(read.csv(inFile, stringsAsFactors=FALSE))
}

simplify <- function(genes, rules){
  srules <- list()
  for (r  in 1:length(rules)) {
    rule <- gsub(" ", "", rules[r])
    rule <- gsub("&", ",",rule, fixed=TRUE)
    rule <- gsub("|", ",",rule, fixed=TRUE)
    rule <- gsub("(", "", rule, fixed=TRUE)
    rule <- gsub(")", "", rule, fixed=TRUE)
    rule <- gsub("!","", rule, fixed=TRUE)
    srules[r] <- list(which(genes %in% strsplit(rule,",")[[1]]))
  }
  return (srules)
}


countDegree <- function(network, name, noOfSims = 500, noOfThreads=8 ) {
  ccNw <- readNetwork(network, name)
  ccGenes <- ccNw[,1]
  ccRules <- ccNw[,2]
  sr <- simplify(ccGenes,ccRules)
  degs <- matrix(0,ncol=6, nrow=length(ccGenes), dimnames=list(ccGenes,list("out","inp","loop","total","Z (out)", "Z (total)")))
  degs[,1] <- sapply(1:length(ccGenes),function(g) { length(sr[[g]]) })
  degs[,2] <- c(tabulate(unlist(sr)), rep(0,length(ccGenes)-max(unlist(sr))))
  degs[,3] <- sapply(1:length(ccGenes), function(g) { if (g %in% sr [[g]]) {return(1)} else {return(0)} } )
  degs[,1] <- degs[,1] - degs[,3]
  degs[,2] <- degs[,2] - degs[,3]
  degs[,4] <-  degs[,1] + degs[,2] + degs[,3]
  degMean <- mean(degs[,1])
  degSd <- sd(degs[,1])
  degs[,5] <- round((degs[,1] - degMean)/degSd,2)
  degMean <- mean(degs[,4])
  degSd <- sd(degs[,4])
  degs[,6] <- round((degs[,4] - degMean)/degSd,2)
  degtable <- degs
  degtable[degtable==0] <- NA
  m <- degs[,4]
  m <- m[m>0]
  m_pl <- displ$new(m)
  est <- estimate_xmin(m_pl)
  if(is.na(est$xmin)) {
    print("failed to set model parameters")
    pvalue <- "failed to set model parameters"
  } else {
    m_pl$setXmin(est)
    bt_pl <- bootstrap_p(m_pl, no_of_sims=noOfSims, threads=noOfThreads)
    pvalue <- bt_pl$p
  }
  return(list(bootstrap=bt_pl, degrees=degs))
}

network <- "PanNet_Model_NoDelays.txt"
Zscore <- countDegree(network,network, noOfSims=100, noOfThreads=4)
print(Zscore)

##scale freeness result
#$p
#[1] 0.4
#Note, that genes/proteins with a score Z(total)>2.5 are defined as hub nodes.

########Figure S1#######
conv2adjmat <- function(sbmlnet, inputcorrected = FALSE){
  "Converts an SBML object to adjacency matrix. 
  If inputcorrected = T, all input edges are set to zero. 
  A vertex v is said to be an input if it is only regulated by itself, 
  meaning the sum of column v in the adjacency matrix is one, with the only non-zero entry 
  being at position [v,v]."
  adjmat <- sapply(sbmlnet$interactions, function(gene) {v <- rep(0,length(sbmlnet$genes)); 
  v[gene$input] <- 1; return(v)})
  if (inputcorrected == TRUE){
    for (d in 1:dim(adjmat)[1]){
      if (adjmat[d,d] >= 1){adjmat[d,d] <- 1}
    }
    for (col in 1:dim(adjmat)[2]){
      if (sum(adjmat[,col]) == 1 & adjmat[col,col] == 1){
        adjmat[col,col] <-0
      }
    }
  }
  return(adjmat)
}

PanNET_adjmat <- conv2adjmat(PanNET)
PanNET_noDelays_adjmat <- conv2adjmat(PanNET_noDelays)

######Figure S1a - Interaction graph based on Boolean functions with node size related to Z-scores
totaldegs_noDelays <- rep(NA, length(PanNET_noDelays$genes))
for (g in 1:length(PanNET_noDelays$genes)){
  totaldegs_noDelays[g] <- sum(PanNET_noDelays_adjmat[g,]) + sum(PanNET_noDelays_adjmat[,g])
}
#hist(totaldegs_noDelays, breaks = seq(1,18,1), main = "Total degree distribution in network without time delays")

totaldegs_noDelays_zscores <- rep(NA, length(PanNET_noDelays$genes))
for (g in 1:length(PanNET_noDelays$genes)){
  totaldegs_noDelays_zscores[g] <- (totaldegs_noDelays[g] - mean(totaldegs_noDelays))/sd(totaldegs_noDelays)
}

g <- graph_from_adjacency_matrix(PanNET_noDelays_adjmat)

#####Coloring according to Figure 3a
#Angiogenesis
redgenes <- c("VEGF", "VEGFR", "ET1", "ARF1", "STAT3", "HIF1") #HIF1 green-red, STAT3 red-blue, ARF1 red-purple
#mTORC signaling
bluegenes <- c("S6K", "TSC", "mTORC1", "mTORC2", "FourEBP1")
#MAPK signaling
purplegenes <- c("FAK1", "RAS", "RAF", "MEK", "ERK", "NEK2")
#PI3K/AKT Signaling
yellowgenes <- c("GLI1", "AMPK", "RSK", "IGF", "IGFR", "IRS", "PP2A", "CIP2A", "AP1", #yellow-purple
                 "PI3K", "PTEN", "AKT", "FOXO3", "FOXO1", "SKP2", "GSK3B")
#Cell cycle
greengenes <- c("eIFs", "REDD", "SETD7", #green-blue
                "RASSF1A", #green-yellow
                "CTNNB1", "PDX1", #green-grey
                "SLUG1", "TP53", "MDM2", "DAXX", "E2F", "CCND1", "RB", "P18", "MEN1", "P21", "MYC", "CCNE1", "P27")
#Cell adhesion
greygenes <- c("RAC1", "CDH1", "RAC1", "IQGAP", "TJ")

#rest, i.e. delay nodes remain white
colorvec <- rep("white", length(PanNET_noDelays$genes))
colorvec[which(PanNET_noDelays$genes %in% redgenes)] <- "#E2A096"
  colorvec[which(PanNET_noDelays$genes %in% bluegenes)] <- "#9BCDD6"
    colorvec[which(PanNET_noDelays$genes %in% purplegenes)] <- "#DAB8FF"
      colorvec[which(PanNET_noDelays$genes %in% yellowgenes)] <- "#F7DD93"
        colorvec[which(PanNET_noDelays$genes %in% greengenes)] <- "#B0DC00"
          colorvec[which(PanNET_noDelays$genes %in% greygenes)] <- "#D6D6D6"
            
          nodelabelsize <- (2*totaldegs_noDelays_zscores+10)*0.05
          l <- readRDS("./layoutcoords_interactiongraph.RDS")
          plot(g, layout=l, vertex.size=2*totaldegs_noDelays_zscores+10, vertex.label.color="black",
               vertex.color=colorvec, edge.color="grey", vertex.label.cex=nodelabelsize,
               edge.arrow.size=0.3)
          
####Figure S1b  
 #Total degree = sum up both row and col for each gene to get a vector          
 totaldegs_noDelays <- rep(NA, length(PanNET_noDelays$genes))
 for (g in 1:length(PanNET_noDelays$genes)){
   totaldegs_noDelays[g] <- sum(PanNET_noDelays_adjmat[g,]) + sum(PanNET_noDelays_adjmat[,g])
 }
# hist(totaldegs_noDelays, breaks = seq(1,18,1), main = "Total degree distribution in network without time delays")
 
 totaldegs_noDelays_zscores <- rep(NA, length(PanNET_noDelays$genes))
 for (g in 1:length(PanNET_noDelays$genes)){
   totaldegs_noDelays_zscores[g] <- (totaldegs_noDelays[g] - mean(totaldegs_noDelays))/sd(totaldegs_noDelays)
 } 
 
 ggplot(as.data.frame(totaldegs_noDelays), aes(x=totaldegs_noDelays)) + 
   geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Total node degree distribution in network without time delays") +
   xlab("Total node degree") + ylab("Frequency") + 
   theme_bw() + scale_y_continuous(breaks=seq(0,10,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,18,1)) + 
   theme(panel.grid.minor = element_blank())
 
 ####Figure S1c          
 #In-degree = sum up cols across the adjmat to get a vector  
 indegs_noDelays <- colSums(PanNET_noDelays_adjmat)
# hist(indegs_noDelays) 
 ggplot(as.data.frame(indegs_noDelays), aes(x=indegs_noDelays)) + 
   geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Node in-degree distribution in network without time delays") +
   xlab("Node in-degree") + ylab("Frequency") + 
   theme_bw() + scale_y_continuous(breaks=seq(0,21,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,11,1)) + 
   theme(panel.grid.minor = element_blank())
 
 ####Figure S1d    
 #Out-degree = sum up rows across the adjmat to get a vector
 outdegs_noDelays <- rowSums(PanNET_noDelays_adjmat)
# hist(outdegs_noDelays)
 
 ggplot(as.data.frame(outdegs_noDelays), aes(x=outdegs_noDelays)) + 
   geom_histogram(binwidth = 1, color="black", fill="grey") + ggtitle("Node out-degree distribution in network without time delays") +
   xlab("Node out-degree") + ylab("Frequency") + 
   theme_bw() + scale_y_continuous(breaks=seq(0,19,1), expand=c(0,0)) + scale_x_continuous(breaks=seq(1,13,1)) + 
   theme(panel.grid.minor = element_blank())
 ############################################################################################################################################
 ######Figure S2 - Stability assessment

 modTestNetworkProperties <- function (network, numRandomNets = 100, testFunction = "testIndegree", 
                                       testFunctionParams = list(), accumulation = c("characteristic", 
                                                                                     "kullback_leibler"), alternative = c("greater", "less"), 
                                       sign.level = 0.05, drawSignificanceLevel = TRUE, klBins, 
                                       klMinVal = 1e-05, linkage = c("uniform", "lattice"), functionGeneration = c("uniform", 
                                                                                                                   "biased"), validationFunction, failureIterations = 10000, 
                                       simplify = FALSE, noIrrelevantGenes = TRUE, d_lattice = 1, 
                                       zeroBias = 0.5, title = "", xlab, xlim, breaks = 30, ...) 
 {
   stopifnot(inherits(network, "BooleanNetwork") || inherits(network, 
                                                             "SymbolicBooleanNetwork"))
   if (is.character(testFunction)) 
     testFunctionName <- testFunction
   else testFunctionName <- ""
   testFunction <- match.fun(testFunction)
   accumulate <- (match.arg(accumulation) == "characteristic")
   origResult <- testFunction(network, accumulate, testFunctionParams)
   numGenes <- length(network$interactions)
   if (inherits(network, "SymbolicBooleanNetwork")) 
     inputGenes <- sapply(network$interactions, function(interaction) length(getInputs(interaction)))
   else inputGenes <- sapply(network$interactions, function(interaction) length(interaction$input))
   if (missing(validationFunction)) 
     validationFunction <- NULL
   randomResults <- lapply(seq_len(numRandomNets), function(i) {
     randomNet <- generateRandomNKNetwork(n = numGenes, k = inputGenes, 
                                          topology = "fixed", linkage = linkage, functionGeneration = functionGeneration, 
                                          validationFunction = validationFunction, failureIterations = failureIterations, 
                                          simplify = simplify, noIrrelevantGenes = noIrrelevantGenes, 
                                          d_lattice = d_lattice, zeroBias = zeroBias)
     randomRes <- testFunction(randomNet, accumulate, testFunctionParams)
     return(randomRes)
   })
   if (accumulate) 
     randomResults <- unlist(randomResults)
   args <- list(...)
   res <- switch(match.arg(accumulation, c("characteristic", 
                                           "kullback_leibler")), characteristic = {
                                             if (missing(xlab)) {
                                               xlab <- switch(testFunctionName, testIndegree = "Gini index of state in-degrees", 
                                                              testAttractorRobustness = "% of identical attractors", 
                                                              testTransitionRobustness = "Normalized Hamming distance", 
                                                              "accumulated results")
                                             }
                                             if (missing(xlim)) {
                                               xlim <- range(c(origResult, randomResults))
                                             }
                                             alternative <- match.arg(alternative, c("greater", "less"))
                                             if (alternative == "greater") pval <- sum(randomResults < 
                                                                                         origResult)/length(randomResults) else pval <- sum(randomResults > 
                                                                                                                                              origResult)/length(randomResults)
                                             if (testFunctionName == "testIndegree" | testFunctionName == 
                                                 "testAttractorRobustness") {
                                               r <- hist(randomResults, xlim = xlim, xlab = xlab, 
                                                         main = title, xaxt = "n", ...)
                                               axis(side = 1, at = seq(xlim[1], xlim[2], length.out = 11))
                                             } else {
                                               r <- hist(randomResults, xlim = xlim, xlab = xlab, 
                                                         main = title, ...)
                                             }
                                             print("orig Result :")
                                             print(origResult)
                                             abline(v = origResult, col = "red")
                                             if (alternative == "greater") text(x = origResult, pos = 2, 
                                                                                y = max(r$counts) * 0.75, labels = paste("> ", round(pval * 
                                                                                                                                       100), "%\nof random results", sep = ""), col = "red", 
                                                                                cex = 0.75) else text(x = origResult, pos = 4, y = max(r$counts) * 
                                                                                                        0.75, labels = paste("< ", round(pval * 100), "%\nof random results", 
                                                                                                                             sep = ""), col = "red", cex = 0.75)
                                             if (drawSignificanceLevel) {
                                               if (alternative == "greater") {
                                                 quant <- quantile(randomResults, 1 - sign.level)
                                                 abline(v = quant, col = "blue")
                                                 text(x = quant, pos = 2, y = max(r$counts) * 
                                                        0.85, labels = paste((1 - sign.level) * 100, 
                                                                             "% quantile", sep = ""), col = "blue", cex = 0.75)
                                               } else {
                                                 quant <- quantile(randomResults, sign.level)
                                                 print("random Result :")
                                                 print(quant)
                                                 abline(v = quant, col = "blue")
                                                 text(x = quant, pos = 4, y = max(r$counts) * 
                                                        0.85, labels = paste(sign.level * 100, "% quantile", 
                                                                             sep = ""), col = "blue", cex = 0.75)
                                               }
                                             }
                                             list(hist = r, pval = 1 - pval, significant = (1 - pval <= 
                                                                                              sign.level))
                                           }, kullback_leibler = {
                                             if (missing(xlab)) xlab <- "Kullback-Leibler distance"
                                             if (missing(klBins)) {
                                               bins <- unique(c(origResult, unlist(randomResults)))
                                               bins <- c(bins, max(bins) + 1)
                                             } else {
                                               bins <- unique(c(origResult, unlist(randomResults)))
                                               if (klBins < length(bins)) bins <- seq(min(bins), 
                                                                                      max(bins), length.out = klBins + 1) else bins <- c(bins, 
                                                                                                                                         max(bins) + 1)
                                             }
                                             vals <- sapply(randomResults, function(results) kullbackLeiblerDistance(origResult, 
                                                                                                                     results, bins = bins, minVal = klMinVal))
                                             r <- hist(vals, xlab = xlab, main = title, breaks = breaks, 
                                                       ...)
                                             list(hist = r)
                                           }, stop("'accumulation' must be one of \"characteristic\",\"kullback_leibler\""))
   return(res)
 }
 set.seed(10000)
 a<-loadNetwork(file = "PanNet_10_Mar_NoDelays.txt")
 Hamming<-modTestNetworkProperties(a, numRandomNets = 1000, 
                                   testFunction = "testTransitionRobustness",
                                   testFunctionParams = list(numSamples=1000),
                                   alternative="less")
 
 #results for the Draft in Figure S5
 #[1] "orig Result :"
 #[1] 0.01490909
 #[1] "random Result :"
 #5% 
 #0.02574242 
 
 ############################################################################################################################################
 ######Figure S4
 attr<-simulateSymbolicModel(PanNET, startStates = 10000)
 plotAttractors(attr, allInOnePlot = T, onColor = "olivedrab3", offColor = "slategrey", title = "Unperturbed PNET", drawLabels = F, drawLegend = F)
 #3 attractors --> 1x single state (98.92%), 2x 8-states (0.72%, 0.37%) --> with 10k startStates
 #3 attractors --> 1x single state (99.16), 2x 8-states (0.52%, 0.33%) --> with 1 MIli startStates
 
 ############################################################################################################################################
 ######Figure S5
 DAXXko<-fixGenes(PanNET, "DAXX", 0)
 attrDAXXko<-simulateSymbolicModel(DAXXko, startStates = 10000)
 plotAttractors(attrDAXXko, allInOnePlot = T, onColor = "olivedrab3", offColor = "slategrey", title ="DAXX loss", drawLegend = F)
 #1x single state (98.8%), 1x 2-states (0.09%), 2x 8-states (0.65%, 0.46%)
 
 ############################################################################################################################################
 ######Figure S6
 TSCko<-fixGenes(PanNET, "TSC", 0)
 attrTSCko<-simulateSymbolicModel(TSCko, startStates = 10000)
 plotAttractors(attrTSCko, allInOnePlot = T, onColor = "olivedrab3", offColor = "slategrey", title = "TSC loss", drawLegend = F)
 #3 attractors --> 1x single state (99.25%), 2x 8-states (0.53%, 0.22%)
 
 ############################################################################################################################################
 ######Figure S7
 MEN1ko<-fixGenes(PanNET, "MEN1", 0)
 attrMEN1ko<-simulateSymbolicModel(MEN1ko, startStates = 10000)
 plotAttractors(attrMEN1ko, onColor = "olivedrab3", offColor = "slategrey", title = "MEN1 loss", drawLegend = F)
 #12 attractors--> 1x single state (1.26%), 1x 2-states (20.51%), 2x 3-states (0.02%, 0.34%), 1x 4-states (19.64%), 4x 6-states (13.12%, 0.65%, 0.08%, 0.06%), 2x 8-states (21.72%, 18.46%), 1x 30-states (4.14%)
 
 ############################################################################################################################################
 #Figure S9
 DAXXMENKO<-fixGenes(DAXXko, "MEN1", 0)
 attrDAXXkoMENko<-simulateSymbolicModel(DAXXMENKO, startStates = 10000)
 plotAttractors(attrDAXXkoMENko, allInOnePlot = F, onColor = "olivedrab3", offColor = "slategrey" , title ="DAXX and MEN1 loss", drawLegend = F)

 ############################################################################################################################################
 #FIGURE S11  
 mTORC1ko<-fixGenes(PanNET, "mTORC1", 0)
 attrmTORC1<-simulateSymbolicModel(mTORC1ko, startStates = 10000)
 plotAttractors(attrmTORC1, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 intervention PanNET untperturbed", drawLegend = F, allInOnePlot = T)
 
 ############################################################################################################################################
 #FIGURE S12  
 mTORC1koDAXXko<-fixGenes(mTORC1ko, "DAXX", 0)
 attrmTORC1DAXXko<-simulateSymbolicModel(mTORC1koDAXXko, startStates = 10000)
 plotAttractors(attrmTORC1DAXXko, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 intervention in DAXX loss", drawLegend = F)
 
 ############################################################################################################################################
 #FIGURE S13
 mTORC1koTSCko<-fixGenes(mTORC1ko, "TSC", 0)
 attrmTORC1TSCko<-simulateSymbolicModel(mTORC1koTSCko, startStates = 10000)
 plotAttractors(attrmTORC1TSCko, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 intervention in TSC loss", drawLegend = F)
 
 ############################################################################################################################################
 #FIGURE S14  
 mTORC1koMEN1ko<-fixGenes(mTORC1ko, "MEN1", 0)
 attrmTORC1MEN1ko<-simulateSymbolicModel(mTORC1koMEN1ko, startStates = 10000)
 plotAttractors(attrmTORC1MEN1ko, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 intervention in MEN-1 loss", drawLegend = F)
 
 ############################################################################################################################################
 #FIGURE S15
 mTORC1ko<-fixGenes(PanNET, "mTORC1", 0)
 mTORC1mTORC2ko<-fixGenes(mTORC1ko, "mTORC2", 0)
 attrmTORC1mTORC2<-simulateSymbolicModel(mTORC1mTORC2ko, startStates = 10000)
 plotAttractors(attrmTORC1mTORC2, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 and mTORC2 intervention PanNET untperturbed", drawLegend = F)
 
 ############################################################################################################################################
 #FIGURE S16
 mTORC1komTORC2DAXXko<-fixGenes(mTORC1mTORC2ko, "DAXX", 0)
 attrmTORC1mTORC2DAXXko<-simulateSymbolicModel(mTORC1komTORC2DAXXko, startStates = 10000)
 plotAttractors(attrmTORC1mTORC2DAXXko, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 and mTORC2 intervention in DAXX loss", drawLegend = F)
 
 ############################################################################################################################################
 #FIGURE S17
 mTORC1komTORC2TSCko<-fixGenes(mTORC1mTORC2ko, "TSC", 0)
 attrmTORC1mTORC2TSCko<-simulateSymbolicModel(mTORC1komTORC2TSCko, startStates = 10000)
 plotAttractors(attrmTORC1mTORC2TSCko, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 and mTORC2 intervention in TSC loss", drawLegend = F)
 
 ############################################################################################################################################
 #FIGURE S18
 mTORC1komTORC2MEN1ko<-fixGenes(mTORC1mTORC2ko, "MEN1", 0)
 attrmTORC1mTORC2MEN1ko<-simulateSymbolicModel(mTORC1komTORC2MEN1ko, startStates = 10000)
 plotAttractors(attrmTORC1mTORC2MEN1ko, onColor = "olivedrab3", offColor = "slategrey", title = "mTORC1 and mTORC2 intervention in MEN-1 loss", drawLegend = F)
 
 ############################################################################################################################################
 ###Analysis to generate Figure S19-S20

 MEN1ko<-fixGenes(PanNET, "MEN1", 0)
 
 #exhaustive
 attrMEN1koex<-simulateSymbolicModel(MEN1ko,  method = "sat.exhaustive")
 plotAttractors(attrMEN1koex, onColor = "olivedrab3", offColor = "slategrey", title = "MEN1 loss EX", drawLegend = F)
 
 #1,000 start states
 attrMEN1ko1k<-simulateSymbolicModel(MEN1ko, startStates = 1000)
 plotAttractors(attrMEN1ko1k, onColor = "olivedrab3", offColor = "slategrey", title = "MEN1 loss 1k", drawLegend = F)
 
 #10,000 start states
 attrMEN1ko10k<-simulateSymbolicModel(MEN1ko, startStates = 10000)
 plotAttractors(attrMEN1ko10k, onColor = "olivedrab3", offColor = "slategrey", title = "MEN1 loss 10k", drawLegend = F)
 
 #100,000 start states
 attrMEN1ko100k<-simulateSymbolicModel(MEN1ko, startStates = 100000)
 plotAttractors(attrMEN1ko100k, onColor = "olivedrab3", offColor = "slategrey", title = "MEN1 loss 100k", drawLegend = F)
 
 #1,000,000 start states
 #might require server running
 attrMEN1ko1M<-simulateSymbolicModel(MEN1ko, startStates = 1000000)
 plotAttractors(attrMEN1ko1M, onColor = "olivedrab3", offColor = "slategrey", title = "MEN1 loss 1M", drawLegend = F)
 

#####Results are summarized in BasinAnalysisMeanSD.csv and BasinAnalysisByPhenMeanSD.csv
#####Please note that the randomly sampling may slightly change the results

############################################################################################################################################
######Figure S19 
 BasinDistribution<-read.csv("BasinAnalysisMeanSD.csv", sep = ",", header = T)
 BasinDistribution$Attractor<-factor(BasinDistribution$Attractor, levels = c("Attractor1", "Attractor2", "Attractor3", "Attractor4", "Attractor5", "Attractor6", "Attractor7", "Attractor8" ,"Attractor9", "Attractor10", "Attractor11", "Attractor12", "Attractor13", "Attractor14" ))
 BasinDistribution$States<-factor(BasinDistribution$States, levels = c("1k", "10k", "100k", "1M"))   
 
 ggplot(data=BasinDistribution, aes(x=Attractor, y=Basin_Mean, fill=States)) +
   geom_bar(stat="identity", position=position_dodge()) +
   geom_errorbar(aes(ymin=Basin_Mean-Basin_SD, ymax=Basin_Mean+Basin_SD), width=0.9,
                 position=position_dodge(0.9)) +
   scale_fill_manual(values=c("#96C217", "#2E7DFF", "#CFA7FF","#F5D586")) +
   theme(axis.text.x = element_text(angle=90)) +
   scale_y_continuous(limits =c(0, 50), breaks = seq(0, 50, by = 2))

 ############################################################################################################################################
 ######Figure S20 
 BasinDistributionPhen<-read.csv("BasinAnalysisByPhenMeanSD.csv", sep = ",", header = T)
 BasinDistributionPhen$Phenotype<-factor(BasinDistributionPhen$Phenotype, levels = c("Angio_G0_Prolif", "Quiescence", "G0_alert", "Angio_G0" ))
 BasinDistributionPhen$States<-factor(BasinDistributionPhen$States, levels = c("1k", "10k", "100k", "1M"))   
 
 ggplot(data=BasinDistributionPhen, aes(x=Phenotype, y=Basin_Mean, fill=States)) +
   geom_bar(stat="identity", position=position_dodge()) +
   geom_errorbar(aes(ymin=Basin_Mean-Basin_SD, ymax=Basin_Mean+Basin_SD), width=0.9,
                 position=position_dodge(0.9)) +
   scale_fill_manual(values=c("#96C217", "#2E7DFF", "#CFA7FF","#F5D586")) +
   theme(axis.text.x = element_text(angle=90)) +
   scale_y_continuous(limits =c(0, 100), breaks = seq(0, 100, by = 20))
 
 ############################################################################################################################################ 
 