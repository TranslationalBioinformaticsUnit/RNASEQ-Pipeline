#######################################################
######       DIFFERENTIAL EXPRESSION AND GSEA    ######
#######################################################

#load files
setwd("C:/Users/alberto.maillo/Desktop/CellLines/")
load("RPKM_cqn.Rda")
load("class.Rda")
load("ensembleID.Rda")
load("DifferentialExpresion.Rda")

#############
#1-Differential Expresion (Separate CellLine)
#############

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
#make contrast (11 contrast)
cont.matrix.TIME <- makeContrasts(
  #Cell-line: WM
  WM_Sus5_Adh5= cWM_suspension_5 - cWM_adhesion_5,
  WM_Sus25_Adh25= cWM_suspension_25 - cWM_adhesion_25,
  WM_Adh5_Adh25= cWM_adhesion_5 - cWM_adhesion_25,
  WM_Sus5_Sus25= cWM_suspension_5 - cWM_suspension_25,
  WMefectGlucose=((cWM_suspension_5 - cWM_adhesion_5) - (cWM_suspension_25 - cWM_adhesion_25)),
  #Cell-line: 501
  c501_Sus5_Adh5= c501_suspension_5 - c501_adhesion_5,
  c501_Sus25_Adh25= c501_suspension_25 - c501_adhesion_25,
  c501_Adh5_Adh25= c501_adhesion_5 - c501_adhesion_25,
  c501_Sus5_Sus25= c501_suspension_5 - c501_suspension_25,
  c501efectGlucose=((c501_suspension_5 - c501_adhesion_5) - (c501_suspension_25 - c501_adhesion_25)),
  adhesionEffectCellLine = ((cWM_adhesion_5 - cWM_adhesion_25)-(c501_adhesion_5 - c501_adhesion_25)),
  levels=design)

fitLimma.TIME<- contrasts.fit(fitLimma, cont.matrix.TIME)
fitLimma.TIME <- eBayes(fitLimma.TIME)

#Join all information
Diff_All<- data.frame(1:dim(RPKM.cqn)[1])
for(i in 1:11)
{#i<-1

  Diff_TIMEgo<-topTable(fitLimma.TIME,coef=i,
                        number=nrow(RPKM.cqn),
                        sort.by="none")
  Diff_TIMEgo<-Diff_TIMEgo[,c("logFC","t","P.Value","adj.P.Val")]
  colnames(Diff_TIMEgo) <- paste(colnames(cont.matrix.TIME)[i],colnames(Diff_TIMEgo),sep="_")
  Diff_All <- cbind(Diff_All,Diff_TIMEgo)
}
Diff_All<-Diff_All[,-1]

save(Diff_All, file="DifferentialExpresion.Rda")


###########
#2-GSEA
###########
load("DifferentialExpresion.Rda")
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")

interestingDiff<-c("Adh5_Sus5", "Adh25_Sus25")
cell<-c("WM", "c501")

#Notes:
#WM with LogFC and pValueCutoff=0.05: OK
#501 with LogFC and pValueCutoff=0.05: No terms -> up PVlaueCutOff=0.1: No terms
#-> down minGSize: No terms; with adjP.Val-> No terms; with P.Val->No terms; with t_value: No terms

#GO Gene Set Enrichment Analysis
#cell
set.seed(1234567)
for(i in 1:2){
  #interestingDiff
  for(j in 1:2){
    ##2.1-Prepare geneList (Ranked List)
    #tValue
    geneList<-Diff_All[, paste(cell[i], interestingDiff[j], "t", sep="_")]
    names(geneList)<-rownames(Diff_All)
    geneList<- sort(geneList, decreasing = TRUE)
    
    ##2.2 Gene Set Enrichment Analysis
    gse<-gseGO(geneList = geneList, OrgDb = org.Hs.eg.db,
               
               
               keyType = 'ENSEMBL',
               ont = "ALL", nPerm = 1000,
               minGSSize = 20, maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose=FALSE)
    gse<-setReadable(gse, OrgDb = org.Hs.eg.db)
    write.table(gse,paste("GSEAtVal_", cell[i],"_", interestingDiff[j],".txt", sep=""),
                sep="\t",quote=F,row.names=F)
  }
}


#############
#2ºOption (GO Over Representation)
#cell
for (i in 1:2){
  #interestingDiff
  for(j in 1:2){
    ##2.1-Prepare geneList
    #adj.P.Val
    geneList<-Diff_All[,paste(cell[i], interestingDiff[j], "adj.P.Val", sep="_")]
    names(geneList)<-rownames(Diff_All)
    #genes high differentiated
    gene<-names(geneList[geneList<0.001])
    
    ##2.2 GO over-representation test
    gse2<-enrichGO(gene = gene,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   keyType = 'ENSEMBL',
                   ont='ALL',
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable=TRUE)
    rownames(ensembleID)<-ensembleID[,1]
    genesIdSymbol<-ensembleID[gene,5]
    #delete NA values
    genesIdSymbol<-genesIdSymbol[genesIdSymbol!=""]
    write.table(genesIdSymbol, paste("GenesOverRepre_", cell[i],"_", interestingDiff[j],".txt", sep=""))
    write.table(gse2,paste("OverRepre_", cell[i],"_", interestingDiff[j],".txt", sep=""),
                sep="\t",quote=F,row.names=F)
  }
}


#########
#PLOTS
#######
i<-1 
j<-1
load("C:/Users/alberto.maillo/Desktop/CellLines/GO/gse.Rda")
geneList<-Diff_All[,paste(cell[i], interestingDiff[j], "logFC", sep="_")]
names(geneList)<-rownames(Diff_All)
geneList<- sort(geneList, decreasing = TRUE)

gseaplot(gse, geneSetID = "GO:0022613")
gseaplot(gse, geneSetID = "GO:0006270")
#FEW GENES(20)
gseaplot(gse, geneSetID = "GO:0001556")
library(ggplot2)
gsearank(gse, geneSetID = "GO:0022613")

#####################################################################################################
#####################################################################################################
############################
####Differential expresion
####difference between cell lines- between adhesion and suspension
setwd("C:/Users/alberto.maillo/Desktop/CellLines/")
load("RPKM_cqn.Rda")
load("class.Rda")
load("ensembleID.Rda")

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

DiffCellLineGrowth<-topTable(fitLimma.TIME,coef=1,number=nrow(RPKM.cqn),sort.by="none")
sum(DiffCellLineGrowth$adj.P.Val<0.05) # 9034 genes

save(DiffCellLineGrowth, file="DiffCellLineGrowth.Rda")

#################
#GO difference between cell lines- between adhesion and suspension
################
library("clusterProfiler")
library("org.Hs.eg.db")
library("DOSE")


##Prepare geneList
#tValue
geneList<-DiffCellLineGrowth[,"t"]
names(geneList)<-rownames(DiffCellLineGrowth)
geneList<- sort(geneList, decreasing = TRUE)

##Gene Set Enrichment Analysis
gse<-gseGO(geneList = geneList, OrgDb = org.Hs.eg.db,
           keyType = 'ENSEMBL',
           ont = 'BP', nPerm = 1000,
           minGSSize = 20, maxGSSize = 500,
           pvalueCutoff = 0.05,
           verbose=FALSE)
gse<-setReadable(gse, OrgDb = org.Hs.eg.db)
write.table(gse,paste("GSEA_CellLinestVal_Adh_Sus.txt", sep=""),
            sep="\t",quote=F,row.names=F)
