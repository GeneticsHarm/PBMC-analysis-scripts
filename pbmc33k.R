##
## Script for clustering and assigning cell types to the 33k PBMC dataset.
##
## Loaded data:
##     ~/Documents/Stage/Matrices/count-matrix-33k.tsv
## Saved data:
##     /data/seurat_33k_merged


library(Seurat)
library(Matrix)

# Load data.
pbmc.data <- read.table("~/Documents/Stage/Matrices/count-matrix-33k.tsv", row.names=1, header=T)
# Remove all genes with zero counts.
pbmc.data <- pbmc.data[apply(pbmc.data, 1, function(x) !all(x==0)),]

# Get a dataframe with gene names corresponding to the ensemle id's in the data.
library(biomaRt)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene.name.table <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                         filters = "ensembl_gene_id",
                         values = rownames(pbmc.data),
                         mart = ensembl)

# Not all ensembl ids are present in the BioMart database, so only set when available.
get.gene.name <- function(ensemble.id) {
  gene.name <- gene.name.table[gene.name.table$ensembl_gene_id == ensemble.id, 1]
  if (identical(gene.name, character(0)) || gene.name == "") {
    gene.name <- ensemble.id
  }
  return(gene.name)
}

rownames(pbmc.data) <- make.names(lapply(rownames(pbmc.data), get.gene.name), unique = T)

pbmc.data.sparse <- Matrix(as.matrix(pbmc.data), sparse = T)

## Load this file to start from this point
#save(pbmc.data.sparse, file = "./data/33k_sparse.Rda")
#load("./data/33k_sparse.Rda")

dim(pbmc.data.sparse)

pbmc33k <- new("seurat", raw.data = pbmc.data.sparse)

#setup setting do.scale and do.center to F - this means that we will NOT scale genes by default (to speed things up)
pbmc33k <- Setup(pbmc33k, min.cells = 3, min.genes = 200, project = "10X_PBMC", do.scale = F, do.center = F, names.field = 2, names.delim = "\\.")

length(pbmc33k@cell.names) # 33118 / 33123
dim(pbmc33k@data)
dim(pbmc33k@data) # 23734 33118

mito.genes <- grep("^MT\\.", rownames(pbmc33k@data), value = T)
pbmc33k@data[mito.genes, ]
percent.mito <- colSums(expm1(as.matrix(pbmc33k@data[mito.genes, ]))) / colSums(expm1(as.matrix(pbmc33k@data)))

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
pbmc33k <- AddMetaData(pbmc33k, percent.mito, "percent.mito")
VlnPlot(pbmc33k, c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(pbmc33k, "nUMI", "nGene")

all.cells <- FetchData(pbmc33k, c("nUMI", "nGene", "percent.mito"))
plot(all.cells$nGene~all.cells$nUMI, ylab = "Number of genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.3), cex=0.2, cex.lab=1.4)
abline(h = 2500, lwd = 2, lty = 2, col = "red")
abline(h = 500, lwd = 2, lty = 2, col = "red")

plot(all.cells$percent.mito~all.cells$nUMI, ylab = "Percentage mitochondrial genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.4), cex=0.5, cex.lab=1.4)
abline(h = 0.05, lwd = 2, lty = 2, col = "red")

#We filter out cells that have unique gene counts over 2,500 and under 500, and > 5% mitochondrial percentage
pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.high = 2500)
dim(pbmc33k@data) # 33090/33118 cells remain
pbmc33k <- SubsetData(pbmc33k, subset.name = "nGene", accept.low = 500)
dim(pbmc33k@data) # 30364 cells remain
pbmc33k <- SubsetData(pbmc33k, subset.name = "percent.mito", accept.high = 0.05)
dim(pbmc33k@data) # 29304 cells remain

#choose gene outliers on mean-variability plot
pbmc33k <- MeanVarPlot(pbmc33k.merged, x.low.cutoff = 0, y.cutoff = 0.8)
length(pbmc33k@var.genes)
pbmc33k@var.genes

#Perform negative-binomial regression on the variable genes, this sets their value in pbmc33k@scale.data, which is used for PCA/clustering
#We only do this on the variable genes to save time, but you can do this genome-wide
#We treat mitochondrial percentage, batch, and nUMI as confounding variables, 

#save(pbmc33k, file = "./data/seurat_33k_before_PCA.Rda")
#load("./data/seurat_33k_before_PCA.Rda")c
pbmc33k <- RegressOut(pbmc33k,latent.vars = c("percent.mito", "orig.ident", "nUMI"), genes.regress = pbmc33k@var.genes)

#Run PCA with the IRLBA package (iteratively computes the top dimensions, dramatic increase in speed since we are throwing away most PCs anyway)
pbmc33k <- PCAFast(pbmc33k, pc.genes = pbmc33k@var.genes, pcs.compute = 40, pcs.print = 30)
#pbmc33K <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)


PCElbowPlot(pbmc33k.merged, num.pc = 40)
PrintPCA(pbmc33k,pcs.print = 1:36)
PCAPlot(pbmc33k.merged, 1, 2, cols.use=rep(rgb(0,0,0,alpha=0.4), 28), pt.size = 0.5, no.legend = T)
PCHeatmap(pbmc33k,pc.use = 1:12,100)
PCHeatmap(pbmc33k,pc.use = 13:24,100)
PCHeatmap(pbmc33k,pc.use = 25:36,100)

pbmc33k@var.genes

#select 25 PCs for downstream analysis
pbmc33k <- RunTSNE(pbmc33k, dims.use = 1:25, do.fast = T)

#save.SNN means that you can easily re-run with different resolution values. Here we run with a few different res v
pbmc33k <- FindClusters(pbmc33k ,pc.use = 1:25, resolution = seq(2,4,0.5), save.SNN = T, do.sparse = T)

#pbmc33k_20=FindClusters(pbmc33k,pc.use = 1:25,resolution = seq(0.6,4,0.1),save.SNN = T,do.sparse = T,k.param = 20)
#save(pbmc33k,file = "~/Projects/datasets/pbmc33k/pbmc33k_nbinom_final.Robj")

#Explore results with a resolution value of 2
#We see nice clusters, but the data is slightly under-clustered (for example CD1C+ and CD141+ DCs are merged together)
pbmc33k <- SetAllIdent(pbmc33k, id = "res.2")
TSNEPlot(pbmc33k,do.label = F)


#In the section below, we explore an overclustering combined with post-hoc merging strategy that can help discover weaker splits in the data

#We can bump the resolution up to call more clusters, but this slightly over-clusters the data
pbmc33k <- SetAllIdent(pbmc33k, id = "res.4")
TSNEPlot(pbmc33k)

#In fact there is no 'perfect' value of resolution, we always either under or over-cluster the data
#This is because of dramatically different cluster sizes, and is known as the 'multi-resolution' problem in graph-based clustering
#One solution is to slightly over-cluster the data, and then perform a post-hoc merging step, where transcriptionally indistinguishable clusters are merged back together
#As a test for merging, we use the Out-of-bag error (OOBE) from a random forest classifier, but you could also set a cutoff for # of differentially expressed genes

pbmc33k <- SetAllIdent(pbmc33k, id = "res.4")

#Build a classification hierarchy that places transcriptionally similar clusters adjacent on a tree
pbmc33k <- BuildClusterTree(pbmc33k, do.reorder = T, reorder.numeric = T)

#calculate the classification error for left/right cells on each branch of the tree
#sort internal nodes based on OOBE. For nodes with high OOBE, we cannot accurately tell the left/right children apart based on random forests, so the clusters may need to be merged
node.scores <- AssessNodes(pbmc33k)
node.scores[order(node.scores$oobe,decreasing = T),] -> node.scores

#choose the first eight splits to merge )
nodes.merge=node.scores[1:8,]
nodes.to.merge=sort(nodes.merge$node)
pbmc33k.merged <- pbmc33k
for (n in nodes.to.merge){
  pbmc33k.merged <- MergeNode(pbmc33k.merged, n)
}

#rebuild the classification hierarchy using all genes for interpretation (you can also try with variable genes as the default)
pbmc33k.merged <- BuildClusterTree(pbmc33k.merged, do.reorder = T, reorder.numeric = T,genes.use = rownames(pbmc33k.merged@data))
TSNEPlot(pbmc33k.merged, do.label = F)

#Rebuild the tree to update the leaf identities (upcoming version will update this automatically)
pbmc33 <- BuildClusterTree(pbmc33k.merged, genes.use = rownames(pbmc33k.merged@data))

PlotClusterTree(pbmc33k.merged)

#color TSNE based on a hierarchical split
ColorTSNESplit(pbmc33k.merged, node = 33)
ColorTSNESplit(pbmc33k.merged, node = 31,color1 = "red",color3="blue")

#Visualize canonical markers
FeaturePlot(pbmc33k.merged, c("nUMI"), cols.use = c("lightgrey","blue"))
#naive/memory
FeaturePlot(pbmc33k.merged, c("CCR7"), cols.use = c("lightgrey","blue"))
FeaturePlot(pbmc33k.merged, c("S100A4"), cols.use = c("lightgrey","blue"))

FeaturePlot(pbmc33k.merged, c("CD34", "CD1C", "FCER1A", "CST3"), cols.use = c("lightgrey","blue"), nCol = 2)
FeaturePlot(pbmc33k.merged, c("PPBP"), cols.use = c("lightgrey","blue"))
FeaturePlot(pbmc33k.merged, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"),cols.use = c("lightgrey","blue"),nCol = 3)

dim(pbmc33k.merged@data)

layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

par(mfrow=c(2,2))
VlnPlot(pbmc33k.merged, c("MS4A1"), size.use = 0.2) 
VlnPlot(pbmc33k.merged, c("CD79A"), size.use = 0.2) 
VlnPlot(pbmc33k.merged, c("CD8A"), size.use = 0.2) 
VlnPlot(pbmc33k.merged, c("CD3E"), size.use = 0.05, size.x.use = .1) 
grep("^CD8", rownames(pbmc33k@data), value = T)

cluster3 <- WhichCells(pbmc33k.merged, ident = 3)
plasma <- WhichCells(pbmc33k.merged, ident = "B_Cells_Plasma")
plasma

FindMarkersNode(pbmc33.merged, node = 31, min.diff.pct = 0.25)

#Dot plot visualization
DotPlot(pbmc33k.merged)
DotPlot(pbmc33k.merged,c("CD3E","CD8B","SELL","GNLY","GZMA","GZMB","GZMH","GZMK","PRF1","NKG7","XCL2","FGFBP2","KLRC1"),cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)
# B cell comparison
DotPlot(pbmc33k.merged,c("CD3E","CD79A", "MS4A1", "JCHAIN", "IGHA1", "IGHA2", "IGLC2"), cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2, ylab="Cell type", xlab="Marker gene")
# Cluster 3 comparison
DotPlot(pbmc33k.merged,c("CD3E","GNLY","FGFBP2","GZMA","GZMK","GZMH","GZMB","STMN1","TUBA1B","FCER1A","TYMS","HMGB2"),cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)
DotPlot(pbmc33k.merged,c("CD4", "HLA.DRA", "HLA.DRB5", "HLA.DRB1","IL3RA", "CLEC4C", "TLR7", "ITGAX", "CD14", "THBD", "CD1C"),cols.use = myPalette(low="lightgrey",high = "blue"),cex.use = 2)
grep("^CD1", rownames(pbmc33k@data), value = T)
#pDC markers IL3RA(CD123), CLEC4C(BDCA-2), TLR7.. Negative markers: ITGAX CD 14

pbmc33k.downsampled <- SubsetData(pbmc33k.merged, max.cells.per.ident = 100)
markers <- FindAllMarkersNode(pbmc33k.downsampled, node = 31, max.cells.per.ident = 100)
HeatmapNode(pbmc33k.downsampled, node = 34, marker.list = markers, max.genes = 10)
FindMarkers(pbmc33k.downsampled, ident.1 = 17,  only.pos = T)

#save(pbmc33k.merged, file = "./data/seurat_33k_merged")
#save(pbmc33k, file = "./data/seurat_33k_final.Rda")
#save(pbmc33k.downsampled, file = "./data/seurat_33k_downsampled.Rda")
#load("./data/seurat_33k_final.Rda")
#load("./data/seurat_33k_merged")

#rename cluster IDs
new_ids=c("Megakaryocytes","Monocytes",3,"NK_1","NK_2","B_Cells_Plasma","Dendritic_Cells","B_Cells_1","B_Cells_2","B_Cells_3","T_Cells_CD4+_1","T_Cells_CD8+_1","T_Cells_CD4+_2","T_Cells_CD8+_2","T_Cells_CD8+_3","T_Cells_CD8+_4")
current.cluster.ids=1:16
pbmc33k.merged@ident <- plyr::mapvalues(pbmc33k.merged@ident, from = current.cluster.ids, to = new_ids)

##Find markers for CD1C vs CD141 DCs, using the negative binomial test
##Can speed up if desired  by setting max.cells.per.ident or min.diff.pct
dc.markers=FindMarkers(pbmc33k.merged,"DC_CD1C+","DC_CD141+",test.use = "negbinom")

markers.3 = FindMarkers(pbmc33k.merged, ident.1 = c(3), min.pct = 0.25)
markers.3[1:30,]
markers.6 = FindMarkers(pbmc33k.merged, ident.1 = c(6), min.pct = 0.25)
markers.6[1:30,]

pbmc33k <- pbmc33k.own
load("./data/seurat_33k_merged")
pbmc33k.satija <- pbmc33k

TSNEPlot(pbmc33k.merged, do.label = F, pt.size = 0.5, label.size = 5)
TSNEPlot(pbmc33k.satija, do.label = T)
TSNEPlot(pbmc33k, do.label = T)
TSNEPlot(pbmc33k, do.label = F , pt.size = 0.5)


