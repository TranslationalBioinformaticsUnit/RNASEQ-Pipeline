###############################################
############Visualization Excel################ 
###############################################

setwd("C:/Users/alberto.maillo/Desktop/RNASeq/CreateExcel/")

#################################
#Results visualization 501
#################################

class$glucose <- as.factor(class$glucose)

class501 <- class[which(class$cell_type=="501"),]
m501 <- RPKM.cqn[,which(class$cell_type=="501")]

orderSamples <- order(class501$type_growth, class501$glucose)

intensity_matrix <- m501 - apply(m501,1,median)
write.table(intensity_matrix[,orderSamples],"intensity_matrix_501.tsv", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)

g1 <- Diff_group1c501
g2 <- Diff_group2c501
g3 <- Diff_group3c501
g4 <- Diff_group4c501
g5 <- Diff_group5c501

colnames(g1) <- paste("A5_S5",colnames(g1),sep=":")
colnames(g2) <- paste("A25_S25",colnames(g2),sep=":")
colnames(g3) <- paste("A25_A5",colnames(g3),sep=":")
colnames(g4) <- paste("S25_S5",colnames(g4),sep=":")
colnames(g5) <- paste("(A5_S5)-(A25_S25)",colnames(g5),sep=":")


#gene_symbol <- ensembleID$hgnc_symbol
#gene_symbol[which(ensembleID$hgnc_symbol=="")] <- "NA"
#gene_description<-ensembleID$description
#gene_description[which(ensembleID$description=="")]<-"NA"

gene_info <- cbind(rownames(g1),gene_symbol,g1[,c(1,4,5)],g2[,c(1,4,5)],g3[,c(1,4,5)],g4[,c(1,4,5)],g5[,c(1,4,5)], apply(m501,1,median), apply(m501,1,min), apply(m501,1,max), gene_description)
colnames(gene_info)[1:2] <- c("ensemble_ID","gene_ID")
colnames(gene_info)[18:21] <- c("median","min","max","description")

write.table(gene_info,"gene_info_501.tsv", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)


metadata <- t(class501[orderSamples,c('type_growth','glucose')])
metadata <- cbind(metadata,c("#D","#D"),c("#adhesion=orange","#25=grey"),c("#suspension=blue","#5=yellow"))
write.table(metadata,"metadata_501.tsv", sep="\t", dec=".", row.names=TRUE, col.names=FALSE, quote=FALSE)

rm(g1, g2, g3, g4, g5, class501, m501, orderSamples, intensity_matrix,gene_info, metadata)

##########################################################################
#################################
#Results visualization WM
#################################

class$glucose <- as.factor(class$glucose)

classWM <- class[which(class$cell_type=="WM"),]
mWM <- RPKM.cqn[,which(class$cell_type=="WM")]

orderSamplesWM <- order(classWM$type_growth, classWM$glucose)

intensity_matrixWM <- mWM - apply(mWM,1,median)
write.table(intensity_matrixWM[,orderSamples],"intensity_matrix_WM.tsv", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)

g1WM <- Diff_group1WM
g2WM <- Diff_group2WM
g3WM <- Diff_group3WM
g4WM <- Diff_group4WM
g5WM <- Diff_group5WM

colnames(g1WM) <- paste("A5_S5",colnames(g1WM),sep=":")
colnames(g2WM) <- paste("A25_S25",colnames(g2WM),sep=":")
colnames(g3WM) <- paste("A25_A5",colnames(g3WM),sep=":")
colnames(g4WM) <- paste("S25_S5",colnames(g4WM),sep=":")
colnames(g5WM) <- paste("(A5_S5)-(A25_S25)",colnames(g5WM),sep=":")

# SigG1 <- rep(0,nrow(g1))
# SigG1[which(g1$`A5_S5:logFC`>=0 & g1$`A5_S5:adj.P.Val`<=0.05)] <- "UP"
# SigG1[which(g1$`A5_S5:logFC`<0 & g1$`A5_S5:adj.P.Val`<=0.05)] <- "DW"
# SigG1[which(g1$`A5_S5:logFC`>=0.58 & g1$`A5_S5:adj.P.Val`<=0.05)] <- "UUP"
# SigG1[which(g1$`A5_S5:logFC`<=-0.58 & g1$`A5_S5:adj.P.Val`<=0.05)] <- "DDW"


#gene_symbol <- ensembleID$hgnc_symbol
#gene_symbol[which(ensembleID$hgnc_symbol=="")] <- "NA"
#gene_description<-ensembleID$description
#gene_description[which(ensembleID$description=="")]<-"NA"

gene_infoWM <- cbind(rownames(g1WM),gene_symbol,g1WM[,c(1,4,5)],g2WM[,c(1,4,5)],g3WM[,c(1,4,5)],g4WM[,c(1,4,5)],g5WM[,c(1,4,5)], apply(mWM,1,median), apply(mWM,1,min), apply(mWM,1,max), gene_description)
colnames(gene_infoWM)[1:2] <- c("ensemble_ID","gene_ID")
colnames(gene_infoWM)[18:21] <- c("median","min","max","description")

write.table(gene_infoWM,"gene_info_WM.tsv", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)


metadataWM <- t(classWM[orderSamplesWM,c('type_growth','glucose')])
metadataWM <- cbind(metadataWM,c("#D","#D"),c("#adhesion=orange","#25=grey"),c("#suspension=blue","#5=yellow"))
write.table(metadataWM,"metadata_WM.tsv", sep="\t", dec=".", row.names=TRUE, col.names=FALSE, quote=FALSE)

rm(g1WM, g2WM, g3WM, g4WM, g5WM, classWM, mWM, orderSamplesWM, intensity_matrixWM,gene_infoWM, metadataWM, gene_symbol, gene_description)
#################################################################################
#################################################################################


library("gplots")
#heatmap
PRUEBA<-DiffCellLineGrowth[order(DiffCellLineGrowth$adj.P.Val),]
heatmap.2(RPKM.cqn[rownames(PRUEBA[1:100,]),], col=bluered(245), key=TRUE, key.xlab ="Expression", scale="row", trace="none", xlab= "Samples", ylab="Genes", cexRow=0.25, cexCol = 0.75, srtCol=45, margins = c(6, 8), main="Growth Effect Cell Lines", distfun=function(x) as.dist(1-cor(t(x))))
#############################################################################
#############################################################################

###################
###CREAR EXCEL
#################

Diff_cellLinesAdh
DiffCellLineGrowth

class$glucose <- as.factor(class$glucose)

orderSamples <- order(class$cell_type, class$type_growth, class$glucose)

intensity_matrix <- RPKM.cqn - apply(RPKM.cqn,1,median)
write.table(intensity_matrix[,orderSamples],"intensity_matrix.tsv", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)

g1 <- Diff_cellLinesAdh
g2 <- DiffCellLineGrowth

colnames(g1) <- paste("WM vs 501 A5_A25",colnames(g1),sep=":")
colnames(g2) <- paste("WM vs 501 (suspension vs adhesion)",colnames(g2),sep=":")

gene_info <- cbind(rownames(g1),gene_symbol,g1[,c(1,4,5)],g2[,c(1,4,5)], apply(RPKM.cqn,1,median), apply(RPKM.cqn,1,min), apply(RPKM.cqn,1,max), gene_description)
colnames(gene_info)[1:2] <- c("ensemble_ID","gene_ID")
colnames(gene_info)[9:12] <- c("median","min","max","description")

write.table(gene_info,"gene_info.tsv", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)

metadata <- t(class[orderSamples,c('cell_type','type_growth','glucose')])
metadata <- cbind(metadata,c("#D","#D","#D"),c("#501=red","#adhesion=orange","#25=grey"),c("#WM=green","#suspension=blue","#5=yellow"))
write.table(metadata,"metadata.tsv", sep="\t", dec=".", row.names=TRUE, col.names=FALSE, quote=FALSE)
