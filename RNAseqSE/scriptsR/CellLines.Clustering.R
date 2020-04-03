##############################
######    CLUSTERING    ######
##############################

#load files
setwd("C:/Users/alberto.maillo/Desktop/CellLines/")
load("RPKM_cqn.Rda")
load("class.Rda")
load("ensembleID.Rda")
load("DifferentialExpresion.Rda")


###############################
#1-Obtain genes adj<0.01
#(only for contrasts Adh5vsSus5, Adh25vsSus25)
###############################
Diff_AllWM <- Diff_All[, grepl("WM", colnames(Diff_All))
                       & grepl("Sus", colnames(Diff_All))
                      & grepl("Adh", colnames(Diff_All))
                      & grepl("_adj", colnames(Diff_All)),]

Diff_All501 <- Diff_All[, grepl("c501", colnames(Diff_All))
                       & grepl("Sus", colnames(Diff_All))
                       & grepl("Adh", colnames(Diff_All))
                       & grepl("_adj", colnames(Diff_All)),]

#############1.1-WM
genesHighDifferentiatedWM<-matrix()
#genes <0.001 
#and be high differentiated in both contrast
for(i in 1:ncol(Diff_AllWM))
{
  genesHighDifferentiatedWM <- cbind(genesHighDifferentiatedWM,t(rownames(Diff_AllWM[Diff_AllWM[,i]<0.001,])))
}
genesHighDifferentiatedWM<-genesHighDifferentiatedWM[-1]
genesHighDifferentiatedWM<-genesHighDifferentiatedWM[!duplicated(genesHighDifferentiatedWM)]
length(genesHighDifferentiatedWM) #8336

############1.2-501
genesHighDifferentiated501<-matrix()
#genes <0.001 
#and be high differentiated in both contrast
for(i in 1:ncol(Diff_All501))
{
  genesHighDifferentiated501 <- cbind(genesHighDifferentiated501,t(rownames(Diff_All501[Diff_All501[,i]<0.001,])))
}
genesHighDifferentiated501<-genesHighDifferentiated501[-1]
genesHighDifferentiated501<-genesHighDifferentiated501[!duplicated(genesHighDifferentiated501)]
length(genesHighDifferentiated501) #7333
######



########################
#2-Obtain number of cluster
########################
set.seed(123456)
############2.1-WM
genesHighDifferentiatedWM.norm<-RPKM.cqn[genesHighDifferentiatedWM,which(class$cell_type=="WM")]
#scaled per gene
genesHighDifferentiatedWM.norm.scaled<-t(scale(t(genesHighDifferentiatedWM.norm)))
#numer of cluster
wssWM <- (nrow((genesHighDifferentiatedWM.norm.scaled))-1)*sum(apply(genesHighDifferentiatedWM.norm.scaled, 2, var))
for (i in 2:10){
  wssWM[i] <- sum(kmeans(genesHighDifferentiatedWM.norm.scaled, centers=i)$withinss)
  
}
plot(1:10, wssWM, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares", main="Genes Cluster WM")

###########2.2-501
genesHighDifferentiated501.norm<-RPKM.cqn[genesHighDifferentiated501,which(class$cell_type=="501")]
#scaled per gene
genesHighDifferentiated501.norm.scaled<-t(scale(t(genesHighDifferentiated501.norm)))
#numer of cluster
wss501 <- (nrow((genesHighDifferentiated501.norm.scaled))-1)*sum(apply(genesHighDifferentiated501.norm.scaled, 2, var))
for (i in 2:10){
  wss501[i] <- sum(kmeans(genesHighDifferentiated501.norm.scaled, centers=i)$withinss)
  
}
plot(1:10, wss501, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares", main="Genes Cluster 501")



########################
#3-Clustering
########################
set.seed(123456)

###########3.1-WM
fitKmeans<-kmeans(genesHighDifferentiatedWM.norm.scaled, 2, n=20)
# append cluster assignment
mydataClusterWM <- data.frame(genesHighDifferentiatedWM.norm.scaled, fitKmeans$cluster)
# get cluster means
aggregate(genesHighDifferentiatedWM.norm.scaled, by=list(fitKmeans$cluster),FUN=mean)
colnames(mydataClusterWM)[dim(mydataClusterWM)[2]]<-"Cluster"
sum(mydataClusterWM$Cluster==1)#4443
sum(mydataClusterWM$Cluster==2)#3893
genesHighDifferentiatedWM.norm.scaled <- genesHighDifferentiatedWM.norm.scaled[order(mydataClusterWM$Cluster),]
COLOR_CLUSTERS_WM <- c(rep("green", 4443), rep("yellow",3893))

###########3.2-501
fitKmeans<-kmeans(genesHighDifferentiated501.norm.scaled, 2, n=20)
# append cluster assignment
mydata2Cluster501 <- data.frame(genesHighDifferentiated501.norm.scaled, fitKmeans$cluster)
# get cluster means
aggregate(genesHighDifferentiated501.norm.scaled, by=list(fitKmeans$cluster),FUN=mean)
colnames(mydata2Cluster501)[dim(mydata2Cluster501)[2]]<-"Cluster"
genesHighDifferentiated2501.norm.scaled <- genesHighDifferentiated501.norm.scaled[order(mydata2Cluster501$Cluster),]
COLOR_CLUSTERS2_501 <- c(rep("green", 3770), rep("yellow",3563))
sum(mydata2Cluster501$Cluster==1)#3770
sum(mydata2Cluster501$Cluster==2)#3563

#4Clusters
fitKmeans<-kmeans(genesHighDifferentiated501.norm.scaled, 4, n=20)
# append cluster assignment
mydata4Cluster501 <- data.frame(genesHighDifferentiated501.norm.scaled, fitKmeans$cluster)
# get cluster means
aggregate(genesHighDifferentiated501.norm.scaled, by=list(fitKmeans$cluster),FUN=mean)
colnames(mydata4Cluster501)[dim(mydata4Cluster501)[2]]<-"Cluster"
genesHighDifferentiated4501.norm.scaled <- genesHighDifferentiated501.norm.scaled[order(mydata4Cluster501$Cluster),]
COLOR_CLUSTERS4_501 <- c(rep("green", 1692), rep("yellow",1494), rep("orange",1871), rep("red",2276))
sum(mydata4Cluster501$Cluster==1)#1692
sum(mydata4Cluster501$Cluster==2)#1494
sum(mydata4Cluster501$Cluster==3)#1871
sum(mydata4Cluster501$Cluster==4)#2276


save(genesHighDifferentiatedWM.norm.scaled, file="WMScaled.Rda")
save(genesHighDifferentiated2501.norm.scaled, file="501Scaled2.Rda")
save(mydataClusterWM, file="mydataClusterWM.Rda")
save(mydata2Cluster501, file="mydata2Cluster501.Rda")
save(genesHighDifferentiated4501.norm.scaled, file="501Scaled4.Rda")
save(mydata4Cluster501, file="mydata4Cluster501.Rda")
#############
#4-Plots
#############
library("gplots")

############4.1-WM
COLOR_CLUSTERS_WM <- c(rep("green", 4443), rep("yellow",3893))
pdf("WMClusterHeatmap.pdf")
heatmap.2(genesHighDifferentiatedWM.norm.scaled,col=bluered(245),dendrogram="column", labRow=NA,
          scale = "row",trace = "none", Rowv = FALSE, margins = c(8, 6),cexCol = 0.75,
          xlab = "Samples", ylab="Genes",main="Clusters WM", RowSideColors=COLOR_CLUSTERS_WM)
legend("left", legend = c("cluster1", "cluster2"),
       col = c("green", "yellow"), 
       lty= 1, lwd = 5, cex=.7)
dev.off()

############4.2-501
COLOR_CLUSTERS2_501 <- c(rep("green", 3770), rep("yellow",3563))
pdf("501Cluster2Heatmap.pdf")
heatmap.2(genesHighDifferentiated2501.norm.scaled,col=bluered(245),dendrogram="column", labRow=NA,
          scale = "row",trace = "none", Rowv = FALSE, margins = c(8, 6),cexCol = 0.75,
          xlab = "Samples", ylab="Genes",main="Clusters 501", RowSideColors=COLOR_CLUSTERS2_501)
legend("left", legend = c("cluster1", "cluster2"),
       col = c("green", "yellow"), 
       lty= 1, lwd = 5, cex=.7)
dev.off()
#4Clusters
COLOR_CLUSTERS4_501 <- c(rep("green", 1692), rep("yellow",1494), rep("orange",1871), rep("red",2276))
pdf("501Cluster4Heatmap.pdf")
heatmap.2(genesHighDifferentiated4501.norm.scaled,col=bluered(245),dendrogram="column", labRow=NA,
          scale = "row",trace = "none", Rowv = FALSE, margins = c(8, 6),cexCol = 0.75,
          xlab = "Samples", ylab="Genes",main="Clusters 501", RowSideColors=COLOR_CLUSTERS4_501)
legend("left", legend = c("cluster1", "cluster2", "cluster3", "cluster4"),
       col = c("green", "yellow", "orange", "red"), 
       lty= 1, lwd = 5, cex=.7)
dev.off()


###################################
#5-GSEA ClusterProfile
#################################
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")
library("enrichplot")

##########5.1-WM
#until nClusters
load("mydataClusterWM.Rda")
for (i in 1:2)
{#i<-1
  clusterGenesWM<- rownames(mydataClusterWM[mydataClusterWM$Cluster==i,])
  #funtion enrichCluster
  enrichClusterWM<- enrichGO(gene=clusterGenesWM, OrgDb= org.Hs.eg.db,
                           keyType='ENSEMBL',
                           ont= 'BP',
                           pAdjustMethod = 'BH',
                           pvalueCutoff= 0.01,
                           qvalueCutoff= 0.01)
  
  enrichClusterWM <- setReadable(enrichClusterWM, OrgDb = org.Hs.eg.db)
  #export to txt files
  genesIdSymbol<-ensembleID[clusterGenesWM,5]
  #delete NA values
  genesIdSymbol<-genesIdSymbol[genesIdSymbol!=""]
  write.table(genesIdSymbol, paste("Clustering/WMFINALgenesClusterList_", i, ".txt", sep=""))
  write.table(enrichClusterWM,paste("Clustering/WMGSEA_clusterEnrich_", i, ".txt", sep=""),
              sep="\t",quote=F,row.names=F)
  #save file
  save(enrichClusterWM, file=paste("WMenrichCluster", i, ".Rda", sep=""))
  message(paste("Done with WM cluster ", i, sep = ""))
}

#Compare Clusters WM
clusterGenesWM1<-rownames(mydataClusterWM[mydataClusterWM$Cluster==1,])
clusterGenesWM2<-rownames(mydataClusterWM[mydataClusterWM$Cluster==2,])
clusterGenesTotal<-list(clusterGenesWM1, clusterGenesWM2)
names(clusterGenesTotal)<-c("cluster1", "cluster2")

ckWM<-compareCluster(geneCluster = clusterGenesTotal, fun = "enrichGO", OrgDb= org.Hs.eg.db,
                   keyType='ENSEMBL',
                   ont= 'BP',
                   pAdjustMethod = 'BH',
                   pvalueCutoff= 0.01,
                   qvalueCutoff= 0.01,
                   readable=TRUE)
write.table(ckWM,"CompareClusterEnrich.txt",
            sep="\t",quote=F,row.names=F)

#plots of WM

for(i in 1:2){
 load(paste("Clustering/WMenrichCluster", i, ".Rda", sep = ""))
  if(dim(enrichClusterWM)[1]>0){
    #BARPLOT
    pdf(paste("Clustering/WMbarplotCluster", i, ".pdf",sep=""), width=11)
    barplot(enrichClusterWM, drop=TRUE, showCategory=20)
    dev.off()
    #DOTPLOT
    pdf(paste("Clustering/WMdotplotCluster", i, ".pdf",sep=""), width=11)
    dotplot(enrichClusterWM, showCategory=20)
    dev.off()
    #EMAPPLOT
    pdf(paste("Clustering/WMemapplotCluster", i, ".pdf",sep=""), width=12)
    emapplot(enrichClusterWM)
    dev.off()
    #CNETPPLOT
    pdf(paste("Clustering/WMcnetplotCluster", i, ".pdf",sep=""), width=12)
    cnetplot(enrichClusterWM, categorySize="pvalue")
    dev.off()
  }
}
##########5.2-501
#until nClusters
load("mydata2Cluster501.Rda")
for (i in 1:2)
{#i<-1
  clusterGenes501<- rownames(mydata2Cluster501[mydata2Cluster501$Cluster==i,])
  #funtion enrichCluster
  enrich2Cluster501<- enrichGO(gene=clusterGenes501, OrgDb= org.Hs.eg.db,
                             keyType='ENSEMBL',
                             ont= 'BP',
                             pAdjustMethod = 'BH',
                             pvalueCutoff= 0.01,
                             qvalueCutoff= 0.01)
  
  enrich2Cluster501 <- setReadable(enrich2Cluster501, OrgDb = org.Hs.eg.db)
  #export to txt files
  rownames(ensembleID)<-ensembleID[,1]
  genesIdSymbol<-ensembleID[clusterGenes501,5]
  #delete NA values
  genesIdSymbol<-genesIdSymbol[genesIdSymbol!=""]
  write.table(genesIdSymbol, paste("Clustering/501FINALgenes2ClusterList_", i, ".txt", sep=""))
  write.table(enrich2Cluster501,paste("Clustering/501GSEA_cluster2Enrich_", i, ".txt", sep=""),
              sep="\t",quote=F,row.names=F)
  #save file
  save(enrich2Cluster501, file=paste("501enrich2Cluster", i, ".Rda", sep=""))
  message(paste("Done with 501 cluster ", i, sep = ""))
}

#Compare 2 Clusters 501mel 
cluster2Genes501mel1<-rownames(mydata2Cluster501[mydata2Cluster501$Cluster==1,])
cluster2Genes501mel2<-rownames(mydata2Cluster501[mydata2Cluster501$Cluster==2,])
clusterGenesTotal<-list(cluster2Genes501mel1, cluster2Genes501mel2)
names(clusterGenesTotal)<-c("cluster1", "cluster2")

ck501mel2<-compareCluster(geneCluster = clusterGenesTotal, fun = "enrichGO", OrgDb= org.Hs.eg.db,
                     keyType='ENSEMBL',
                     ont= 'BP',
                     pAdjustMethod = 'BH',
                     pvalueCutoff= 0.01,
                     qvalueCutoff= 0.01,
                     readable=TRUE)
write.table(ck501mel2,"Compare2ClusterEnrich501.txt",
            sep="\t",quote=F,row.names=F)


#until nClusters
load("mydata4Cluster501.Rda")
for (i in 1:4)
{#i<-1
  clusterGenes501<- rownames(mydata4Cluster501[mydata4Cluster501$Cluster==i,])
  #funtion enrichCluster
  enrich4Cluster501<- enrichGO(gene=clusterGenes501, OrgDb= org.Hs.eg.db,
                               keyType='ENSEMBL',
                               ont= 'BP',
                               pAdjustMethod = 'BH',
                               pvalueCutoff= 0.01,
                               qvalueCutoff= 0.01)
  
  enrich4Cluster501 <- setReadable(enrich4Cluster501, OrgDb = org.Hs.eg.db)
  #export to txt files
  rownames(ensembleID)<-ensembleID[,1]
  genesIdSymbol<-ensembleID[clusterGenes501,5]
  #delete NA values
  genesIdSymbol<-genesIdSymbol[genesIdSymbol!=""]
  write.table(genesIdSymbol, paste("Clustering/501FINALgenes4ClusterList_", i, ".txt", sep=""))
  write.table(ensembleID[clusterGenes501,5], paste("Clustering/501genes4ClusterList_", i, ".txt", sep=""))
  write.table(enrich4Cluster501,paste("Clustering/501GSEA_cluster4Enrich_", i, ".txt", sep=""),
              sep="\t",quote=F,row.names=F)
  #save file
  save(enrich4Cluster501, file=paste("501enrich4Cluster", i, ".Rda", sep=""))
  message(paste("Done with 501 cluster ", i, sep = ""))
}

#Compare 4 Clusters 501mel 
cluster4Genes501mel1<-rownames(mydata4Cluster501[mydata4Cluster501$Cluster==1,])
cluster4Genes501mel2<-rownames(mydata4Cluster501[mydata4Cluster501$Cluster==2,])
cluster4Genes501mel3<-rownames(mydata4Cluster501[mydata4Cluster501$Cluster==3,])
cluster4Genes501mel4<-rownames(mydata4Cluster501[mydata4Cluster501$Cluster==4,])
clusterGenesTotal<-list(cluster4Genes501mel1, cluster4Genes501mel2, cluster4Genes501mel3, cluster4Genes501mel4)
names(clusterGenesTotal)<-c("cluster1", "cluster2", "cluster3", "cluster4")

ck501mel4<-compareCluster(geneCluster = clusterGenesTotal, fun = "enrichGO", OrgDb= org.Hs.eg.db,
                          keyType='ENSEMBL',
                          ont= 'BP',
                          pAdjustMethod = 'BH',
                          pvalueCutoff= 0.01,
                          qvalueCutoff= 0.01,
                          readable=TRUE)
write.table(ck501mel4,"Compare4ClusterEnrich501.txt",
            sep="\t",quote=F,row.names=F)



#PLOTS
#nCluster=2
for(i in 1:2){
  load(paste("Clustering/501enrich2Cluster", i, ".Rda", sep = ""))
  if(dim(enrich2Cluster501)[1]>0){
    #BARPLOT
    pdf(paste("Clustering/501barplot2Cluster", i, ".pdf",sep=""), width=11)
    barplot(enrich2Cluster501, drop=TRUE, showCategory=20)
    dev.off()
    #DOTPLOT
    pdf(paste("Clustering/501dotplot2Cluster", i, ".pdf",sep=""), width=11)
    dotplot(enrich2Cluster501, showCategory=20)
    dev.off()
    #EMAPPLOT
    pdf(paste("Clustering/501emapplot2Cluster", i, ".pdf",sep=""), width=12)
    emapplot(enrich2Cluster501)
    dev.off()
    #CNETPPLOT
    pdf(paste("Clustering/501cnetplot2Cluster", i, ".pdf",sep=""), width=12)
    cnetplot(enrich2Cluster501, categorySize="pvalue")
    dev.off()
  }
}

#nCluster=4
for(i in 1:4){
  load(paste("Clustering/501enrich4Cluster", i, ".Rda", sep = ""))
  if(dim(enrich4Cluster501)[1]>0){
    #BARPLOT
    pdf(paste("Clustering/501barplot4Cluster", i, ".pdf",sep=""), width=11)
    barplot(enrich4Cluster501, drop=TRUE, showCategory=20)
    dev.off()
    #DOTPLOT
    pdf(paste("Clustering/501dotplot4Cluster", i, ".pdf",sep=""), width=11)
    dotplot(enrich4Cluster501, showCategory=20)
    dev.off()
    #EMAPPLOT
    pdf(paste("Clustering/501emapplot4Cluster", i, ".pdf",sep=""), width=12)
    emapplot(enrich4Cluster501)
    dev.off()
    #CNETPPLOT
    pdf(paste("Clustering/501cnetplot4Cluster", i, ".pdf",sep=""), width=12)
    cnetplot(enrich4Cluster501, categorySize="pvalue")
    dev.off()
  }
}

#####################################################################################################
#####################################################################################################
############################
####Clustering
####clustering between cell lines- between adhesion and suspension
setwd("C:/Users/alberto.maillo/Desktop/CellLines/")
load("DiffCellLineGrowth.Rda")
##############
#1-Take genes high differentiated (CellLines)
##############
genesHighDiffCellLines<-rownames(DiffCellLineGrowth[DiffCellLineGrowth[,"adj.P.Val"]<0.001,])
length(genesHighDiffCellLines) #6201

########################
#2-Obtain number of cluster (CellLines)
########################
set.seed(123456)
genesHighDiffCellLines.norm<-RPKM.cqn[genesHighDiffCellLines,]
#scaled per gene
genesHighDiffCellLines.norm.scaled<-t(scale(t(genesHighDiffCellLines.norm)))
#numer of cluster
wss <- (nrow((genesHighDiffCellLines.norm.scaled))-1)*sum(apply(genesHighDiffCellLines.norm.scaled, 2, var))
for (i in 2:10){
  wss[i] <- sum(kmeans(genesHighDiffCellLines.norm.scaled, centers=i)$withinss)
  
}
plot(1:10, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares", main="Genes Cluster CellLines")

nClusters<-6


########################
#3-Clustering (CellLines)
########################
set.seed(123456)
fitKmeans<-kmeans(genesHighDiffCellLines.norm.scaled, nClusters, n=20)
# append cluster assignment
mydataClusterCellLines <- data.frame(genesHighDiffCellLines.norm.scaled, fitKmeans$cluster)
# get cluster means
aggregate(genesHighDiffCellLines.norm.scaled, by=list(fitKmeans$cluster),FUN=mean)
colnames(mydataClusterCellLines)[dim(mydataClusterCellLines)[2]]<-"Cluster"
sum(mydataClusterCellLines$Cluster==1)#1510
sum(mydataClusterCellLines$Cluster==2)#1737
sum(mydataClusterCellLines$Cluster==3)#1401
sum(mydataClusterCellLines$Cluster==4)#895
sum(mydataClusterCellLines$Cluster==5)#658
genesHighDiffCellLines.norm.scaled <- genesHighDiffCellLines.norm.scaled[order(mydataClusterCellLines$Cluster),]

save(genesHighDiffCellLines.norm.scaled, file="genesHighDiffCellLines.norm.scaled6.Rda")
#############
#4-Plots (CellLines)
#############
library("gplots")
COLOR_CLUSTERS_CELLLINES<- c(rep("green", 1021), rep("yellow",587), rep("blue", 1006), rep("red",832), rep("orange",1421), rep("black", 1334))
heatmap.2(genesHighDiffCellLines.norm.scaled, col=bluered(245), dendrogram="column", labRow=NA,
          scale = "row",trace = "none", Rowv = FALSE, margins = c(8, 6),cexCol = 0.75,
          xlab = "Samples", ylab="Genes",main="Clusters Cell Lines", RowSideColors=COLOR_CLUSTERS_CELLLINES)

legend("left", legend = c("cluster1", "cluster2", "cluster3", "cluster4","cluster5", "cluster6"),
       col = c("green", "yellow", "blue", "red", "orange", "black"), 
       lty= 1, lwd = 5, cex=.7)

#################################
#5-GSEA ClusterProfile (CellLines)
#################################
for (i in 1:nClusters)
{#i<-1
  clusterGenesCellLines<- rownames(mydataClusterCellLines[mydataClusterCellLines$Cluster==i,])
  #funtion enrichCluster
  enrichClusterCellLines<- enrichGO(gene=clusterGenesCellLines, OrgDb= org.Hs.eg.db,
                              keyType='ENSEMBL',
                              ont= 'BP',
                              pAdjustMethod = 'BH',
                              pvalueCutoff= 0.01,
                              qvalueCutoff= 0.01)
  
  enrichClusterCellLines <- setReadable(enrichClusterCellLines, OrgDb = org.Hs.eg.db)
  #export to txt files
  rownames(ensembleID)<-ensembleID[,1]
  write.table(ensembleID[clusterGenesCellLines,5], paste("Clustering/CellLinesgenesClusterList_", i, "eliminar.txt", sep=""))
  write.table(enrichClusterCellLines,paste("Clustering/CellLinesGSEA_clusterEnrich_", i, "eliminar.txt", sep=""),
              sep="\t",quote=F,row.names=F)
  #save file
  save(enrichClusterCellLines, file=paste("CellLinesenrichCluster", i, "eliminar.Rda", sep=""))
  message(paste("Done with CellLines cluster ", i, sep = ""))
}
