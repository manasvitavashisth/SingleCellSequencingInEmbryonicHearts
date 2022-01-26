library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsignif)

setwd("/Applications/UPenn/SingleCellSequencingHearts/E5")
#setting the working directory

pbmc.data=Read10X(data.dir="/Applications/UPenn/SingleCellSequencingHearts/E5/DMSO/filtered_feature_bc_matrix") #reading the data into R

pbmc.data1=Read10X(data.dir="/Applications/UPenn/SingleCellSequencingHearts/E5/Blebbistatin/filtered_feature_bc_matrix") #reading the data into R
pbmc.data2=Read10X(data.dir="/Applications/UPenn/SingleCellSequencingHearts/E5/BlebbWO/filtered_feature_bc_matrix") #reading the data into R
pbmc.data3=Read10X(data.dir="/Applications/UPenn/SingleCellSequencingHearts/E5/H2O2/filtered_feature_bc_matrix") #reading the data into R

pbmc <- CreateSeuratObject(counts = pbmc.data,project = "DMSO") #creating Seurat object
pbmc1 <- CreateSeuratObject(counts = pbmc.data1,project = "Blebb") #creating Seurat object
pbmc2 <- CreateSeuratObject(counts = pbmc.data2,project = "BlebbWO") #creating Seurat object
pbmc3 <- CreateSeuratObject(counts = pbmc.data3,project = "H2O2") #creating Seurat object

pbmc.combined <- merge(pbmc, y = c(pbmc1,pbmc2,pbmc3), add.cell.ids = c("DMSO", "BLEBB","BLEBBWO","H2O2"), project = "E5hearts")
pbmc=pbmc.combined
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 12)
Idents(pbmc)=pbmc$orig.ident
blebb_markers <- FindMarkers(pbmc, ident.1 = "Blebb", ident.2 = "DMSO")
#change the type of test to test out all of the sets see code below
#FindMarkers(pbmc, ident.1 = "Blebb", ident.2 = "DMSO", test.use = "MAST")

#change the SAMPLE type to test out all of the samples see code below
#FindMarkers(pbmc, ident.1 = "BlebbWO", ident.2 = "DMSO", test.use = "MAST")

genes=c("SMC2", "CENPF", "CENPE", "CDKN1A","TOP2A" )
#change the gene names to your genes of interest
a=blebb_markers[rownames(blebb_markers) %in% genes,]
plot(a$p_val,a$avg_log2FC,xlab="p-value",ylab = "avg_log2FC")
text(a$p_val,a$avg_log2FC,rownames(a))


pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc, features = all.genes)

cluster1.markers <- FindMarkers(pbmc, ident.1 = "BlebbWO",ident.2="DMSO",test.use = "MAST")
a=cluster1.markers(cluster1.markers$genenames=="COL1A1")

cluster1.markers <- FindMarkers(pbmc, ident.1 = "H2O2",ident.2="DMSO",test.use="bimod")
a=cluster1.markers(cluster1.markers$genenames=="COL1A1") 

Idents(pbmc)=pbmc$orig.ident

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mat=as.matrix(pbmc.combined@assays$RNA@counts) 
mat=as.matrix(pbmc[["RNA"]]@data)

pbmc.combined <- merge(pbmc, y = c(pbmc2,pbmc3), add.cell.ids = c("DMSO", "BLEBBWO","H2O2"), project = "E5hearts")
pbmc=pbmc.combined
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 12)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

gene_of_interest=c("LMNA","COL1A1","BRCA1","BRCA2","XRCC5","XRCC4","LMNB1","CDKN1")
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = gene_of_interest, repel = TRUE)
plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:9, cells = 100, balanced = TRUE)

DimPlot(pbmc,dims=c(1,5))

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:15)

JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.1)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims=1:15)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

DimPlot(pbmc, reduction = "umap",label = TRUE)

DotPlot(pbmc,features = c())

Idents(pbmc)=pbmc$orig.ident

Idents(pbmc)=pbmc$seurat_clusters

Idents(pbmc)=pbmc$cell_label
# find all markers of cluster 2
cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)

write.table(pbmc.markers, file = "0_6clusters_allmarkers",  sep = "\t")


A=pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
write.table(A, file = "0_6clusters_10markers",  sep = "\t")

cluster7.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25, logfc.threshold = 0.25)

head(cluster7.markers, n = 5)

write.table(cluster7.markers, file = "5fromOclusters_allmarkers",  sep = "\t")

cluster1.markers <- FindMarkers(pbmc, ident.1 = "BlebbWO",ident.2="DMSO",test.use = "MAST")
a=cluster1.markers(cluster1.markers$genenames=="COL1A1")

cluster1.markers <- FindMarkers(pbmc, ident.1 = "H2O2",ident.2="DMSO",test.use="bimod")
a=cluster1.markers(cluster1.markers$genenames=="COL1A1") 





cluster1.markers <- FindMarkers(pbmc, ident.1 = "Epithelial_H2O2",ident.2="Epithelial_DMSO",test.use="bimod")
a=cluster1.markers(cluster1.markers$genenames=="COL1A1")


cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster1.markers, file = "cluster1_allmarkers",  sep = "\t")

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster0.markers, file = "cluster0_allmarkers",  sep = "\t")
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster2.markers, file = "cluster2_allmarkers",  sep = "\t")
cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster4.markers, file = "cluster4_allmarkers",  sep = "\t")

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster5.markers, file = "cluster5_allmarkers",  sep = "\t")

cluster6.markers <- FindMarkers(pbmc, ident.1 = "Macrophages", min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster6.markers, file = "cluster6_allmarkers",  sep = "\t")

FeaturePlot(pbmc, features = c("COL1A1","FAP","COL5A1","COL1A2","COL3A1","POSTN","THBS2","VCAN","ITGA11","COL6A3"),order = TRUE)


FeaturePlot(pbmc, features = c("LMNA","COL1A1","COL1A2","COL3A1","MMP2","ACTA2","COL4A1","BRCA1","BRCA2","XRCC5","XRCC4","LMNB1","FOXM1","CENPF","SMC2","CDKN1A"),order = TRUE)

FeaturePlot(pbmc,features=c("COL1A1","COL1A2","COL4A1","COL4A2","LMNB1","FOXM1","ACTA2","TAGLN","POSTN","CENPF","TGFB2","VCAN"),order = TRUE,label = TRUE)

FeaturePlot(pbmc, features=c("TEK","PECAM1","CDH5","POSTN","COL1A1","COL3A1","TNNT2","TNNC1","HBE","HBZ"),label = TRUE,order = TRUE)

FeaturePlot(pbmc, features=c("TGFB2","WNT11"))

FeaturePlot(pbmc, features=c("DCN","CDH5","TNNC1","HBZ"))

FeaturePlot(pbmc, features=c("CRABP-I","LMNA","COL1A1","COL1A2","RARA","RARB","RARG","RXRA","CRABP2"),order=TRUE)

iden=rep(c("Epithelial","Erythrocytes","Valve/Chondrocyte/Mesenchymal","Cardiomyocytes","Endothelial","DMSO_Epithelial","Macrophages"),times=a)

cell_label1=cbind(pbmc$orig.ident,pbmc$seurat_clusters)
cell_label1=as.data.frame(cell_label1)
cell_label=cell_label[order(cell_label$V2),]

write.table(label, file = "cell_label",  sep = "\t")

write.table(cell_label, file = "cell_label1",  sep = "\t")

cell_label=read.table("cell_label",header = TRUE,sep="\t",stringsAsFactors = FALSE,na.strings=c("", "NA"))
cell_label=cbind(cell_label,iden)

VlnPlot(pbmc,features="MSX1")

table(cell_label$V1,cell_label$iden)
table(cell_label1$V2)

a=table(cell_label$iden)
x=table(cell_label$V1,cell_label$iden)
barplot(a)

barplot(table(cell_label$iden),xlab="Cell identity",ylab="Number of cells",col=plasma(7))

barplot(x,xlab="Cell identity",ylab="Number of cells",axes = TRUE,axisnames = TRUE,col=plasma(4),legend.text = c("Blebb","BlebbWO","DMSO","H2O2"),angle = 45,ann =TRUE)

names.arg = c("Cardiomyocytes","DMSO_Epithelial","Endothelial","Epithelial","Erythrocytes","Macrophages","Valve/Chondrocyte/Mesenchymal")

barplot(cell_label$iden~cell_label$V1,data = cell_label,xlab="Sample",ylab="Number of cells", col=plasma(7))

barplot(table(cell_label$V1,cell_label$iden),horiz = TRUE)
x=as.data.frame(x)

cell_label=as.data.frame(cell_label)
boxplot(cell_label$V1~cell_label$iden,xlab="Treatment",ylab="Average LMNB1 expression",col=plasma(7))

label=merge(cell_label1,cell_label,by=0,sort=FALSE)

pbmc[["cell_label"]]=label$iden

pbmc[["cell_label"]]=cell_label$iden

VlnPlot(pbmc,features = c("COL1A1","COL3A1","COL1A2"))
FeaturePlot(pbmc,features = c("COL1A1","COL3A1","COL1A2"),label = TRUE,order = TRUE)

VlnPlot(pbmc,features = c("LECT1","CHODL","SOX9","FRZB","VCAN","THBS2"))
VlnPlot(pbmc,features = c("LECT1","CHODL","SOX9","FRZB"))

VlnPlot(pbmc,features = c("HBE","HBBR","HBBA","HBAD","UROD","TSPO2"))

VlnPlot(pbmc,features = c("TNNT2","TNNC1","TNNI1","ACTA1","ACTN2","NKX2-5"))

VlnPlot(pbmc,features = c("TNNT2","TNNC1","TNNI1","ACTA1","ACTN2","NKX2-5"))

VlnPlot(pbmc,features = c("POSTN","PECAM1","KDR","FN1","ENG","COL4A1"))

VlnPlot(pbmc,features = c("TNNT2","TNNC1","TNNI1","ACTA1","ACTN2","NKX2-5"))

VlnPlot(pbmc,features = c("CRABP-I","LMNA","CRABP2","COL1A1","RARG","MAB21L1"))

VlnPlot(pbmc,features = c("CRABP-I","MAB21L1","OSR1","ADAMTS9","GSC","HOXA3"))

VlnPlot(pbmc,features = c("CSF1R","CD74","LCP1","LY86","GSTA3","TIMD4"))

VlnPlot(pbmc,features = c("CDKN1A","LMNA","COL1A1","COL1A2","FOXM1","LMNB1"))

DotPlot(pbmc,features = c("HBE","HBBR","HBBA","HBAD","UROD","TSPO2","TNNT2","TNNC1","TNNI1","ACTA1","ACTN2","NKX2-5","POSTN","PECAM1","KDR","FN1","ENG","COL4A1","CRABP-I","MAB21L1","OSR1","ADAMTS9","GSC","HOXA3","CSF1R","CD74","LCP1","LY86","GSTA3","TIMD4","LECT1","CHODL","SOX9","FRZB","VCAN","THBS2","COL1A1","COL3A1","COL1A2","LMNB1","FOXM1","BRCA1","BRCA2"))+RotatedAxis()

VlnPlot(pbmc,features = "percent.mt")

VlnPlot(pbmc,features = "LMNA")+geom_signif(comparisons = list(c("DMSO","BlebbWO"),c("DMSO","H2O2")),map_signif_level = TRUE,textsize = 5)+ylim(0,6)

VlnPlot(pbmc,features = "SMC2")+geom_signif(comparisons = list(c("DMSO","H2O2")),map_signif_level = TRUE,textsize = 5)+ylim(0,5)

RidgePlot(pbmc,features = "LMNA")

