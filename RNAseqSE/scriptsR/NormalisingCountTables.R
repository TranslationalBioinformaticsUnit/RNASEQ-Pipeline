#source("https://bioconductor.org/biocLite.R")
#biocLite("cqn")
#biocLite("stringi")
#biocLite("biomaRt")
#biocLite("EDASeq")
#biocLite("bit")
#biocLite("edgeR")
#biocLite("BiocUpgrade")
setwd("C:/Users/alberto.maillo/Desktop/CellLines/CountTablesGenes/")
#setwd("C:/Users/alberto.maillo/Desktop/kaustRNAseq/CountTables/")

###########################################################
#######-------UNIFIED COUNT TABLES---------################
###########################################################

#take all the tables files from the directory
count_tables<-dir()

sample_name<-unlist(strsplit(count_tables[1], "_nonempty"))
count_table_unified <- read.table(count_tables[1])
rownames_count_table_unified <- count_table_unified[,1]
count_table_unified <- data.frame(count_table_unified[,-1])
names(count_table_unified) <- sample_name[1]



#import count_tables (unified)
for (i in 2:length(count_tables)){
  sample_name<-unlist(strsplit(count_tables[i], "_nonempty"))
  table<-read.table(count_tables[i])
  #delete the first column(names)
  table<-table[,-1]
  count_table_unified<-cbind(count_table_unified, table)
  colnames(count_table_unified)[i] <- sample_name[1]
}


#paste rownames
rownames(count_table_unified) <- rownames_count_table_unified

#delete innecessary variables
rm(table)
rm(count_tables)
rm(i)
rm(sample_name)
rm(rownames_count_table_unified)


#delete last four files(information of ambiguous...)
count_table_unified<-count_table_unified[1:21504,]
#Read metadata
class <- read.table("C:/Users/alberto.maillo/Desktop/CellLines/metadata.csv", sep=";", header=TRUE)
class$sample_code <- paste(class$FGC.Pool.ID,class$CRG.Sample.ID,class$barcode.seq, sep="_")

#colname = SampleName (Ex. WM-A-5.1)
count_table_unified <- count_table_unified[, paste(0,class$sample_code,sep="")]
colnames(count_table_unified)<-class$Sample.Name
tcount_table_unified <- t(count_table_unified)
all(rownames(tcount_table_unified %in% class$Sample.Name))


###########################################################
#######-------MDS(before normalising)---------#############
###########################################################


#SAMPLES IN ROWNAMES
d<- dist(tcount_table_unified)
fit<- cmdscale(d, eig=TRUE, k=2)

#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

#---------------CELL TYPE(BEFORE CQN)-----------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(Before CQN) by cell type", type="n")
text(x, y, labels = row.names(tcount_table_unified), cex=0.7, col=as.numeric(class$cell_type))
legend("topright",levels(class$cell_type), fill=1:2, cex=0.7)


#--------------GROWTH (BEFORE CQN)--------------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(Before CQN) by growth", type="n")
text(x, y, labels = row.names(tcount_table_unified), cex=0.7, col=as.numeric(class$type_growth))
legend("topright",levels(class$type_growth), fill=1:2, cex=0.7)

#------GLUCOSE CONCENTRATION (BEFORE CQN)-------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(Before CQN) by glucose concentration", type="n")
text(x, y, labels = row.names(tcount_table_unified), cex=0.7, col=as.numeric(as.factor(class$glucose)))
legend("topright",levels(as.factor(class$glucose)), fill=1:2, cex=0.7)


###########################################################
#######----------NORMALISING(CQN)-------------#############
###########################################################

#----FILTERS BEFORE CQN-------------------------------------------
#---1ºFilter:Delete the genes that have not been mapped to any sample
#[18701 x 24] genes
count_table_unified_genes_mapped<-count_table_unified[which(rowSums(as.matrix(count_table_unified))!=0),]
dim(count_table_unified_genes_mapped)
#---2ºFilter:CPM
#------Keep genes that have 1 read per million in at least three samples
library("edgeR")
cpm<-cpm(count_table_unified_genes_mapped)
#[13429x24] genes
count_table_unified_genes_mapped<-count_table_unified_genes_mapped[which((rowSums(cpm>=1))>2), ]
dim(count_table_unified_genes_mapped)

rm(tcount_table_unified, cpm)


#----GET ENSEMBL IDS OF THE GENES----------------
library("biomaRt")
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#--------ensembl-----------------
ensembleID <- getBM(attributes=c('ensembl_gene_id','percentage_gene_gc_content' ,'start_position', 'end_position', 'hgnc_symbol', 'description'), 
                    filters= c('ensembl_gene_id'), values= rownames(count_table_unified_genes_mapped), mart= mart)
#Delete ("ENSG00000276911","ENSG00000278937","ENSG00000283374")-> no EnsemblID database
count_table_unified_genes_mapped<-count_table_unified_genes_mapped[-which(rownames(count_table_unified_genes_mapped) %in% c("ENSG00000276911","ENSG00000278937","ENSG00000283374")),]


#order
ensembleID<-ensembleID[order(ensembleID$ensembl_gene_id),]
count_table_unified_genes_mapped<-count_table_unified_genes_mapped[order(rownames(count_table_unified_genes_mapped)),]
genes_length<-ensembleID$end_position-ensembleID$start_position

#---------CQN ENSEMBL-------------
library("cqn")
count_table_unified.cqn<-cqn(count_table_unified_genes_mapped, lengths = genes_length ,
                             x =ensembleID$percentage_gene_gc_content, 
                             sizeFactors = colSums(as.matrix(count_table_unified_genes_mapped)) , verbose=TRUE)
#----------NORMALIZED VALUES-------
#RPKM(reads per kilobase of exon model per million reads)
RPKM.cqn<- count_table_unified.cqn$y + count_table_unified.cqn$offset

#Save files
setwd("C:/Users/alberto.maillo/Desktop/CellLines/")
save(RPKM.cqn, file="RPKM_cqn.Rda")
save(class, file="class.Rda")
save(ensembleID, file="ensembleID.Rda")
###########################################################
#######-------MDS(AFTER normalising)---------##############
###########################################################

tRPKM.cqn<-t(RPKM.cqn)

#sample information in rows
dnorm<- dist(tRPKM.cqn)
fit<- cmdscale(dnorm, eig=TRUE, k=2)

#plot solution
x<- fit$points[,1]
y<- fit$points[,2]

#---------------CELL TYPE(AFTER CQN)------------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(AFTER CQN) by cell type", type="n")
text(x, y, labels = row.names(tRPKM.cqn), cex=0.7, col=as.numeric(class$cell_type))
legend("topright",levels(class$cell_type), fill=1:2, cex=0.7)
#--------------GROWTH (AFTER CQN)--------------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(AFTER CQN) by growth", type="n")
text(x, y, labels = row.names(tRPKM.cqn), cex=0.7, col=as.numeric(class$type_growth))
legend("topright",levels(class$type_growth), fill=1:2, cex=0.7)

#-------------GLUCOSE CONCENTRATION (AFTER CQN)-------------
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS(AFTER CQN) by glucose concentration", type="n")
text(x, y, labels = row.names(tRPKM.cqn), cex=0.7, col=as.numeric(as.factor(class$glucose)))
legend("topright",levels(as.factor(class$glucose)), fill=1:2, cex=0.7)

#########################
######BOXPLOT############
#########################
boxplot(log2(count_table_unified_genes_mapped), range=0 , main= "Boxplot before CQN", las=2)
boxplot(RPKM.cqn, main="Boxplot after CQN", las=2)

##########################################
####deleted genes per read per million####
##########################################
x1 <- 1:50
y1 <- vector()
for (i in 1:50){
  tt <- apply(cpm,1,FUN=function(x){return(sum(x>i))})
  y1[i] <- sum(tt<4)
}
par(mar=c(5,7,4,2))
plot(x1,y1,type="b", las=1, ylab="DELETED GENES\n", xlab="MINIMUN CPM PER SAMPLE")

#delete innecessary variables
rm(x1)
rm(y1)
rm(tt)


####################################
#####Plot genes tomados de Irene####
####################################
#After CQN
#AXL<-ensembleID[which(ensembleID[,5]=='AXL'),]
#ATF4<-ensembleID[which(ensembleID[,5]=='ATF4'),]
#MITF<-ensembleID[which(ensembleID[,5]=='MITF'),]
##NTRK2 (ENSG00000148053)-> (don't have 1 cpm in at least three samples)
#TRKB<-ensembleID[which(ensembleID[,5]=='NTRK2'),]
#GLUT1<-ensembleID[which(ensembleID[,5]=='SLC2A1'),]
group <- paste(class$cell_type, class$type_growth, class$glucose, sep="_")

names<-c("AXL", "ATF4", "MITF", "SLC2A1")

for (gen in 1:length(names)) {
  pdf(paste(names[gen],".pdf"))
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    plot(RPKM.cqn[ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))],]~as.factor(group), xlab="Samples", ylab="Normalized Counts", cex.axis=0.5, yex.axis=2, main= cbind(names[gen], ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))]))
    stripchart(RPKM.cqn[ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))],]~as.factor(group), 
               vertical = TRUE, method = "jitter",
               pch = 16, col = ifelse(class$cell_type=="501", "red", "green"),
               add = TRUE)

    plot(RPKM.cqn[ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))],class$cell_type=="501"]~as.factor(group[class$cell_type=="501"]), xlab = "501", ylab="Normalized Counts", cex.axis=0.5, yex.axis=2, main= paste("501", names[gen], sep="-") ,names=c("Adhesion_25","Adhesion_5","Suspension_25","Suspension_5"))
    stripchart(RPKM.cqn[ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))],class$cell_type=="501"]~as.factor(group[class$cell_type=="501"]),
               vertical = TRUE, method = "jitter",
               pch = 21, col ="black", bg= "green",
               add = TRUE)

    plot(RPKM.cqn[ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))],class$cell_type=="WM"]~as.factor(group[class$cell_type=="WM"]), xlab = "WM", ylab="Normalized Counts", cex.axis=0.5, yex.axis=2,main= paste("WM", names[gen], sep="-"), names=c("Adhesion_25","Adhesion_5","Suspension_25","Suspension_5"))
    stripchart(RPKM.cqn[ensembleID$ensembl_gene_id[which((ensembleID[,5]==names[gen]))],class$cell_type=="WM"]~as.factor(group[class$cell_type=="WM"]),
               vertical = TRUE, method = "jitter",
               pch = 21, col ="black", bg= "red",
               add = TRUE)

  dev.off()
}