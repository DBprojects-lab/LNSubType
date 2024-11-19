library(ggplot2)
library(RColorBrewer)

LN_TPM=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/LN_TPM_Final_Symbol.txt",header=T,row.names=1,sep="\t")
LN_Info=read.table("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/LN_SubType_Info.txt",header=T,row.names=1,sep="\t")
LN_TPM=LN_TPM[,rownames(LN_Info)]

##############################################################################
####################################     All gene   ##########################
##############################################################################

logTPM=log2(LN_TPM+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(LN_Info))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=LN_Info$GroupByGene,LNStatus=LN_Info$LNGroup)
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=GroupByGene,shape=LNStatus))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))+
scale_color_manual(values=GroupByGene_Colors)
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/PCA/PCA_AllGene.pdf",width=6,height=4)
	print(g)
dev.off()

##############################################################################
######################### cell type specific genes   #########################
##############################################################################
CellTypeMarkers=read.table("D:/06_CRC/LNscRNAseq/CellTypeIdentification/LN.CellType.Marker.txt",header=T,row.names=1,sep="\t")
CellTypeMarkers=CellTypeMarkers[CellTypeMarkers$avg_log2FC>1,]
CellTypeMarkers=CellTypeMarkers[CellTypeMarkers$gene%in%rownames(LN_TPM),]
table(CellTypeMarkers$cluster)

logTPM=log2(LN_TPM[unique(CellTypeMarkers$gene),]+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(LN_Info))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=LN_Info$GroupByGene,LNStatus=LN_Info$LNGroup)
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=GroupByGene,shape=LNStatus))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))+
scale_color_manual(values=GroupByGene_Colors)
pdf("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/PCA/PCA_AllMarkerGene.pdf",width=6,height=4)
	print(g)
dev.off()





for(c in unique(CellTypeMarkers$cluster)){
CellTypeGene=CellTypeMarkers[CellTypeMarkers$cluster==c,"gene"]
CellTypeGeneExpr=LN_TPM[CellTypeGene,]
logTPM=log2(CellTypeGeneExpr+1)
logTPM=logTPM[rowMeans(logTPM)>0,]
logTPM.pca <- prcomp(t(logTPM))
#head(logTPM.pca$x)
all(rownames(logTPM.pca$x)==rownames(LN_Info))
df_pca=data.frame(PC1=logTPM.pca$x[,1],PC2=logTPM.pca$x[,2],GroupByGene=LN_Info$GroupByGene,LNStatus=LN_Info$LNGroup)
GroupByGene_Colors=colorRampPalette(brewer.pal(7, "Set2"))(length(unique(LN_Info$GroupByGene)))
names(GroupByGene_Colors)=sort(unique(LN_Info$GroupByGene))
g=ggplot(df_pca,aes(x=PC1,y=PC2,color=GroupByGene,shape=LNStatus))+ 
geom_point(size=3)+
theme_bw()+
scale_shape_manual(values=c(19, 8))+
scale_color_manual(values=GroupByGene_Colors)

pdf(paste0("D:/06_CRC/bulkRNAseq/BeiJingSample/SubTypeIdentification/PCA/PCA_",c,".pdf",sep=""),width=6,height=4)
	print(g)
dev.off()
}


