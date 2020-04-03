#############################################################
####DIFFERENTIAL EXPRESSION WITH NORMALIZED DATA Limma######
#############################################################
setwd("C:/Users/alberto.maillo/Desktop/CellLines/")
load("RPKM_cqn.Rda")
load("class.Rda")

set.seed(12345678)
library(limma)
#different conditions(cell_type, type_growth, glucose)
condition<-paste(class$cell_type, class$type_growth, class$glucose, sep = "_")
condition<-paste("c", condition, sep="")
#create design
design<-model.matrix(~0 + condition)
#check colnames
colnames(design) <- gsub("condition","",colnames(design))
fitLimma <- lmFit(RPKM.cqn,design)
#make contrast
cont.matrix.TIME <- makeContrasts(
  #Cell-line: WM
  WMgrowth5= cWM_suspension_5 - cWM_adhesion_5,
  WMgrowth25= cWM_suspension_25 - cWM_adhesion_25,
  WMglucoseAd= cWM_adhesion_5 - cWM_adhesion_25,
  WMglucoseSus= cWM_suspension_5 - cWM_suspension_25,
  WMefectGlucose=((cWM_suspension_5 - cWM_adhesion_5) - (cWM_suspension_25 - cWM_adhesion_25)),
  #Cell-line: 501
  c501growth5= c501_suspension_5 - c501_adhesion_5,
  c501growth25= c501_suspension_25 - c501_adhesion_25,
  c501glucoseAd= c501_adhesion_5 - c501_adhesion_25,
  c501glucoseSus= c501_suspension_5 - c501_suspension_25,
  c501efectGlucose=((c501_suspension_5 - c501_adhesion_5) - (c501_suspension_25 - c501_adhesion_25)),
  adhesionEffectCellLine = ((cWM_adhesion_5 - cWM_adhesion_25)-(c501_adhesion_5 - c501_adhesion_25)),
  levels=design)

fitLimma.TIME<- contrasts.fit(fitLimma, cont.matrix.TIME)
fitLimma.TIME <- eBayes(fitLimma.TIME)

#######################
###WM
#adhesion5 vs suspension5
Diff_group1WM<-topTable(fitLimma.TIME,coef=1,number=nrow(RPKM.cqn),sort.by="none")
#adhesion25 vs suspension25
Diff_group2WM<-topTable(fitLimma.TIME,coef=2,number=nrow(RPKM.cqn),sort.by="none")
#adhesion25 vs adhesion5
Diff_group3WM<-topTable(fitLimma.TIME,coef=3,number=nrow(RPKM.cqn),sort.by="none")
#suspension25 vs suspension5
Diff_group4WM<-topTable(fitLimma.TIME,coef=4,number=nrow(RPKM.cqn),sort.by="none")
#adhesion5-suspension5 vs adhesion25-suspension25
Diff_group5WM<-topTable(fitLimma.TIME,coef=5,number=nrow(RPKM.cqn),sort.by="none")

#ADJUSTED P-VALUE <0.05
#Diff_group1WM<-Diff_group1WM[order(Diff_group1WM$adj.P.Val),]
sum(Diff_group1WM$adj.P.Val<0.05) # 10003 genes

#Diff_group2WM<-Diff_group2WM[order(Diff_group2WM$adj.P.Val),]
sum(Diff_group2WM$adj.P.Val<0.05) # 9636 genes

#Diff_group3WM<-Diff_group3WM[order(Diff_group3WM$adj.P.Val),]
sum(Diff_group3WM$adj.P.Val<0.05) # 240 genes

#Diff_group4WM<-Diff_group4WM[order(Diff_group4WM$adj.P.Val),]
sum(Diff_group4WM$adj.P.Val<0.05) # 1320 genes

#Diff_group5WM<-Diff_group5WM[order(Diff_group5WM$adj.P.Val),]
sum(Diff_group5WM$adj.P.Val<0.05) # 366 genes

#######################
###501
#adhesion5 vs suspension5
Diff_group1c501<-topTable(fitLimma.TIME,coef=6,number=nrow(RPKM.cqn),sort.by="none")
#adhesion25 vs suspension25
Diff_group2c501<-topTable(fitLimma.TIME,coef=7,number=nrow(RPKM.cqn),sort.by="none")
#adhesion25 vs adhesion5
Diff_group3c501<-topTable(fitLimma.TIME,coef=8,number=nrow(RPKM.cqn),sort.by="none")
#suspension25 vs suspension5
Diff_group4c501<-topTable(fitLimma.TIME,coef=9,number=nrow(RPKM.cqn),sort.by="none")
#adhesion5-suspension5 vs adhesion25-suspension25
Diff_group5c501<-topTable(fitLimma.TIME,coef=10,number=nrow(RPKM.cqn),sort.by="none")

#ADJUSTED P-VALUE <0.05
#Diff_group1c501<-Diff_group1c501[order(Diff_group1c501$adj.P.Val),]
sum(Diff_group1c501$adj.P.Val<0.05) # 9063 genes

#Diff_group2c501<-Diff_group2c501[order(Diff_group2c501$adj.P.Val),]
sum(Diff_group2c501$adj.P.Val<0.05) # 8692 genes

#Diff_group3c501<-Diff_group3c501[order(Diff_group3c501$adj.P.Val),]
sum(Diff_group3c501$adj.P.Val<0.05) # 4700 genes

#Diff_group4c501<-Diff_group4c501[order(Diff_group4c501$adj.P.Val),]
sum(Diff_group4c501$adj.P.Val<0.05) # 46 genes

#Diff_group5c501<-Diff_group5c501[order(Diff_group5c501$adj.P.Val),]
sum(Diff_group5c501$adj.P.Val<0.05) # 2658 genes

############################
###between different celllines
#adhesion25 vs adhesion5 (between cell lines)
Diff_cellLinesAdh<-topTable(fitLimma.TIME,coef=11,number=nrow(RPKM.cqn),sort.by="none")
Diff_cellLinesAdh<-cellLinesAdh[order(cellLinesAdh$adj.P.Val),]
sum(Diff_cellLinesAdh$adj.P.Val<0.05) # 1705 genes



#SHARED DIFFERENTIAL GENES BETWEEN CELL_LINES
#Group1: adhesion5 vs suspension5
#7111 genes
sum(rownames(Diff_group1WM[which(Diff_group1WM$adj.P.Val<0.05),]) %in% rownames(Diff_group1c501[which(Diff_group1c501$adj.P.Val<0.05),]))
#Group2: adhesion25 vs suspension25
#6716 genes
sum(rownames(Diff_group2WM[which(Diff_group2WM$adj.P.Val<0.05),]) %in% rownames(Diff_group2c501[which(Diff_group2c501$adj.P.Val<0.05),]))
#Group3: adhesion5 vs adhesion25
#164 genes
sum(rownames(Diff_group3WM[which(Diff_group3WM$adj.P.Val<0.05),]) %in% rownames(Diff_group3c501[which(Diff_group3c501$adj.P.Val<0.05),]))
#Group4: suspension5 vs suspension25 
#11 genes
sum(rownames(Diff_group4WM[which(Diff_group4WM$adj.P.Val<0.05),]) %in% rownames(Diff_group4c501[which(Diff_group4c501$adj.P.Val<0.05),]))
#Group5: Group1 vs Group2 
#121 genes
sum(rownames(Diff_group5WM[which(Diff_group5WM$adj.P.Val<0.05),]) %in% rownames(Diff_group5c501[which(Diff_group5c501$adj.P.Val<0.05),]))

######################################################################################################
####################
#####heatmap########
####################
library("gplots")

#WM
heatmap.2(RPKM.cqn[rownames(Diff_group1WM[1:100,]),which((class$glucose=="5") & (class$cell_type=="WM"))], col=bluered(245), 
          key=TRUE, key.xlab ="Expression", scale="row", trace="none", 
          xlab= "Samples", ylab="Genes", cexRow=0.50, cexCol = 0.75, srtCol=45, 
          margins = c(6, 8), main="WM adh5 vs sus5", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group2WM[1:100,]),which((class$glucose=="25") & (class$cell_type=="WM"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="WM adh25 vs sus25", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group3WM[1:100,]),which((class$type_growth=="adhesion") & (class$cell_type=="WM"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="WM adh25 vs adh5", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group4WM[1:100,]),which((class$type_growth=="suspension") & (class$cell_type=="WM"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="WM sus25 vs sus5", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group5WM[1:100,]),which(class$cell_type=="WM")], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="WM Effect Glucose", distfun=function(x) as.dist(1-cor(t(x))))

#heatmap in general all the information
genesDE1WM<-rownames(Diff_group1WM[which(Diff_group1WM$adj.P.Val<0.05),])
genesDE2WM<-rownames(Diff_group2WM[which(Diff_group2WM$adj.P.Val<0.05),])
genesDE3WM<-rownames(Diff_group3WM[which(Diff_group3WM$adj.P.Val<0.05),])
genesDE4WM<-rownames(Diff_group4WM[which(Diff_group4WM$adj.P.Val<0.05),])
genesDE5WM<-rownames(Diff_group5WM[which(Diff_group5WM$adj.P.Val<0.05),])
genesDEWM<-unique(c(genesDE1WM, genesDE2WM, genesDE3WM, genesDE4WM, genesDE5WM))
heatmap.2(RPKM.cqn[genesDEWM, which(class$cell_type=="WM")], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", labRow = FALSE, xlab= "Samples", ylab="Genes", cexCol = 0.75, srtCol=45, margins = c(6, 4), main="WM Heatmap DE", distfun=function(x) as.dist(1-cor(t(x))))

rm(genesDE1WM, genesDE2WM, genesDE3WM, genesDE4WM, genesDE5WM, genesDEWM)


#501
heatmap.2(RPKM.cqn[rownames(Diff_group1c501[1:100,]),which((class$glucose=="5") & (class$cell_type=="501"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.25, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="501 adh5 vs sus5", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group2c501[1:100,]),which((class$glucose=="25") & (class$cell_type=="501"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="501 adh25 vs sus25", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group3c501[1:100,]),which((class$type_growth=="adhesion") & (class$cell_type=="501"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="501 adh5 vs adh25", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group4c501[1:46,]),which((class$type_growth=="suspension") & (class$cell_type=="501"))], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="501 sus5 vs sus25", distfun=function(x) as.dist(1-cor(t(x))))
heatmap.2(RPKM.cqn[rownames(Diff_group5c501[1:100,]),which(class$cell_type=="501")], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.5, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="501 Effect Glucose", distfun=function(x) as.dist(1-cor(t(x))))

#heatmap in general all the information
genesDE1c501<-rownames(Diff_group1c501[which(Diff_group1c501$adj.P.Val<0.05),])
genesDE2c501<-rownames(Diff_group2c501[which(Diff_group2c501$adj.P.Val<0.05),])
genesDE3c501<-rownames(Diff_group3c501[which(Diff_group3c501$adj.P.Val<0.05),])
genesDE4c501<-rownames(Diff_group4c501[which(Diff_group4c501$adj.P.Val<0.05),])
genesDE5c501<-rownames(Diff_group5c501[which(Diff_group5c501$adj.P.Val<0.05),])
genesDEc501<-unique(c(genesDE1c501, genesDE2c501, genesDE3c501, genesDE4c501, genesDE5c501))
heatmap.2(RPKM.cqn[genesDEc501, which(class$cell_type=="501")], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", labRow = FALSE, xlab= "Samples", ylab="Genes", cexCol = 0.75, srtCol=45, margins = c(6, 4), main="501 Heatmap DE", distfun=function(x) as.dist(1-cor(t(x))))


rm(genesDE1c501, genesDE2c501, genesDE3c501, genesDE4c501, genesDE5c501, genesDEc501)


#Heatmap group6-> group3 between cell lines
Diff_cellLinesAdh<-Diff_cellLinesAdh[order(Diff_cellLinesAdh$adj.P.Val),]
heatmap.2(RPKM.cqn[rownames(Diff_cellLinesAdh[1:100,]),which(class$type_growth=="adhesion")], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.25, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="Adhesion Effect Cell Lines", distfun=function(x) as.dist(1-cor(t(x))))


################################################
###----VennDiagram all the different groups---##
################################################
library(VennDiagram)
data_vennWM <- list(g1=rownames(Diff_group1WM[Diff_group1WM$adj.P.Val<0.05,]), g2=rownames(Diff_group2WM[Diff_group2WM$adj.P.Val<0.05,]), g3=rownames(Diff_group3WM[Diff_group3WM$adj.P.Val<0.05,]), g4=rownames(Diff_group4WM[Diff_group4WM$adj.P.Val<0.05,]), g5=rownames(Diff_group5WM[Diff_group5WM$adj.P.Val<0.05,]))
venn.diagram(data_vennWM, filename= 'VennDiagramWM.png',
             output=TRUE,
             fill= c('orange', 'purple', 'green', 'blue', 'red'),
             category = c('Group1', 'Group2', 'Group3','Group4', 'Group5'),
             main="WM cell line",
             cat.col = c('orange', 'purple', 'green', 'blue', 'red')
)
data_venn501 <- list(g1=rownames(Diff_group1c501[Diff_group1c501$adj.P.Val<0.05,]), g2=rownames(Diff_group2c501[Diff_group2c501$adj.P.Val<0.05,]), g3=rownames(Diff_group3c501[Diff_group3c501$adj.P.Val<0.05,]), g4=rownames(Diff_group4c501[Diff_group4c501$adj.P.Val<0.05,]), g5=rownames(Diff_group5c501[Diff_group5c501$adj.P.Val<0.05,]))
venn.diagram(data_venn501, filename= 'VennDiagram501.png',
             output=TRUE,
             fill=c('orange', 'purple', 'green', 'blue', 'red'),
             category = c('Group1', 'Group2', 'Group3','Group4', 'Group5'),
             main="501 cell line",
             cat.col = c('orange', 'purple', 'green', 'blue', 'red')
)
#################################
#################################
#########Clustering##############
#################################

##############
###### WM
##############
#Partitioning
##Determine number of clusters
wssWM <- (nrow((RPKM.cqn))-1)*sum(apply(RPKM.cqn[, which(class$cell_type=="WM")], 2, var))
for (i in 2:10) wssWM[i] <- sum(kmeans(RPKM.cqn[, which(class$cell_type=="WM")], centers=i)$withinss)

plot(1:10, wssWM, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares", main="WM Genes Cluster")

#K-means cluster
fitKmeansWM <- kmeans(RPKM.cqn[, which(class$cell_type=="WM")], 4) # 4 cluster solution
# get cluster means
aggregate(RPKM.cqn[,which(class$cell_type=="WM")],by=list(fitKmeansWM$cluster),FUN=mean)
# append cluster assignment
mydataClusterWM <- data.frame(RPKM.cqn[,which(class$cell_type=="WM")], fitKmeansWM$cluster) 

##################
##WM
d<- dist(RPKM.cqn[, which(class$cell_type=="WM")])
fit<- cmdscale(d, eig=TRUE, k=2)

#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

#---------------WM genes cluster-----------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="WM genes", type="n")
text(x, y, cex=0.7, col=as.numeric(fitKmeansWM$cluster))
legend("topright",levels(as.factor(fitKmeansWM$cluster)), fill=1:4, cex=0.7)

##############
###### 501
##############
#Partitioning
##Determine number of clusters
wss501 <- (nrow((RPKM.cqn))-1)*sum(apply(RPKM.cqn[, which(class$cell_type=="501")], 2, var))
for (i in 2:10) wss501[i] <- sum(kmeans(RPKM.cqn[, which(class$cell_type=="501")], centers=i)$withinss)

plot(1:10, wss501, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares", main="501 Genes Cluster")

#K-means cluster
fitKmeans501 <- kmeans(RPKM.cqn[, which(class$cell_type=="501")], 4) # 4 cluster solution
# get cluster means
aggregate(RPKM.cqn[,which(class$cell_type=="501")],by=list(fitKmeans501$cluster),FUN=mean)
# append cluster assignment
mydataCluster501 <- data.frame(RPKM.cqn[,which(class$cell_type=="501")], fitKmeans501$cluster) 


######
#plot cluster
d<- dist(RPKM.cqn[, which(class$cell_type=="501")])
fit<- cmdscale(d, eig=TRUE, k=2)

#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

#---------------501 genes cluster-----------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="501 genes", type="n")
text(x, y, cex=0.7, col=as.numeric(fitKmeans501$cluster))
legend("topright",levels(as.factor(fitKmeans501$cluster)), fill=1:4, cex=0.7)

#####################################################################################################
############################
####Differential expresion
####difference between cell lines- between adhesion and suspension

set.seed(12345678)
library(limma)
#different conditions(cell_type, type_growth, glucose)
condition<-paste(class$cell_type, class$type_growth, sep = "_")
condition<-paste("c", condition, sep="")
#create design
design<-model.matrix(~0 + condition)
#check colnames
colnames(design) <- gsub("condition","",colnames(design))
fitLimma <- lmFit(RPKM.cqn,design)
#make contrast
cont.matrix.TIME <- makeContrasts(
  CellLineAdhvsSusp=((cWM_suspension - cWM_adhesion) - (c501_suspension - c501_adhesion)),
  levels=design)

fitLimma.TIME<- contrasts.fit(fitLimma, cont.matrix.TIME)
fitLimma.TIME <- eBayes(fitLimma.TIME)

#######################
#heatmap

DiffCellLineGrowth<-topTable(fitLimma.TIME,coef=1,number=nrow(RPKM.cqn),sort.by="none")
sum(DiffCellLineGrowth$adj.P.Val<0.05) # 9034 genes

