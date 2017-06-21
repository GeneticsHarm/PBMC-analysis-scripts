library(Seurat)
library("Matrix.utils")

pbmc8k.data <- Read10X(data.dir = "~/Documents/Stage/Data/pbmc8k/outs/filtered_gene_bc_matrices/GRCh38/")
colnames(pbmc8k.data) <- paste0(colnames(pbmc8k.data), "_8K")
pbmc6k.data <- Read10X(data.dir = "~/Documents/Stage/Data/pbmc6k/outs/filtered_gene_bc_matrices_mex/hg19/")
colnames(pbmc6k.data) <- paste0(colnames(pbmc6k.data), "_6K")
pbmc.merged.data <- merge.Matrix(pbmc8k.data, pbmc6k.data, by.x=rownames(pbmc8k.data), by.y=rownames(pbmc6k.data), all.x=FALSE, all.y=FALSE)

zero.values.6k <- apply(pbmc6k.data, 1, function(x) sum(x==0))
zero.values.6k.frac <- zero.values.6k / ncol(pbmc6k.data)
zero.values.8k <- apply(pbmc8k.data, 1, function(x) sum(x==0))
zero.values.8k.frac <- zero.values.8k / ncol(pbmc8k.data)

hist(zero.values.8k.frac, breaks = 100, xlim = c(0.8, 1), ylim = c(0,5000))
hist(zero.values.6k.frac, breaks = 100, xlim = c(0.8, 1), ylim = c(0,5000))

pbmc.merged.test <- merge.Matrix(pbmc8k.data, pbmc6k.data, by.x=rownames(pbmc8k.data), by.y=rownames(pbmc6k.data), all.x=T, all.y=T)
pbmc.merged.test[1:10,8000:8010]
pbmc.merged.drop <- drop0(pbmc.merged.test)
pbmc.merged.drop[1:10,1:10]

pbmc.merged <- new("seurat", raw.data = pbmc.merged.data)
pbmc.merged <- Setup(pbmc.merged, min.cells = 3, min.genes = 200, project = "10X_PBMC", do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")

mito.genes <- grep("^MT-", rownames(pbmc.merged@data), value = T)
percent.mito <- colSums(expm1(pbmc.merged@data[mito.genes, ])) / colSums(expm1(pbmc.merged@data))
pbmc.merged <- AddMetaData(pbmc.merged, percent.mito, "percent.mito")

all.cells <- FetchData(pbmc.merged, c("nUMI", "nGene", "percent.mito"))
plot(all.cells$nGene~all.cells$nUMI, ylab = "Number of genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.3), cex=0.2, cex.lab=1.4)
abline(h = 2500, lwd = 2, lty = 2, col = "red")
abline(h = 500, lwd = 2, lty = 2, col = "red")

plot(all.cells$percent.mito~all.cells$nUMI, ylab = "Percentage mitochondrial genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.4), cex=0.5, cex.lab=1.4)
abline(h = 0.05, lwd = 2, lty = 2, col = "red")

pbmc.merged <- SubsetData(pbmc.merged, subset.name = "nGene", accept.high = 2500)
dim(pbmc.merged@data) 
pbmc.merged <- SubsetData(pbmc.merged, subset.name = "nGene", accept.low = 500)
dim(pbmc.merged@data) 
pbmc.merged <- SubsetData(pbmc.merged, subset.name = "percent.mito", accept.high = 0.05)
dim(pbmc.merged@data) 

pbmc.merged <- MeanVarPlot(pbmc.merged, x.low.cutoff = 0, y.cutoff = 0.8)
length(pbmc.merged@var.genes)

pbmc.merged <- RegressOut(pbmc.merged, latent.vars = c("percent.mito", "nUMI"), genes.regress = pbmc.merged@var.genes)
pbmc.merged <- PCAFast(pbmc.merged, pc.genes = pbmc.merged@var.genes, pcs.compute = 40, pcs.print = 30)


PCElbowPlot(pbmc.merged, num.pc = 40)
pbmc.merged <- RunTSNE(pbmc.merged, dims.use = 1:20, do.fast = T)

pbmc.merged <- FindClusters(pbmc.merged, pc.use = 1:20, resolution = 4, save.SNN = T, do.sparse = T)

TSNEPlot(pbmc.merged, do.label = F)
pbmc.merged <- BuildClusterTree(pbmc.merged, do.reorder = T, reorder.numeric = T)

is.6k <- vector(mode = "logical", length = ncol(pbmc.merged@data))
is.6k[grep(".+6K$", colnames(pbmc.merged@data))] <- T
is.6k

FeaturePlot(pbmc.merged, c("nUMI"), cols.use = c("lightgrey","blue"))
FeaturePlot(pbmc.merged, c("nGene"), cols.use = c("lightgrey","blue"))
plot(pbmc.merged@tsne.rot$tSNE_1, pbmc.merged@tsne.rot$tSNE_2, cex = 0.2, col = ifelse(is.6k ,'red','black'))

FeaturePlot(pbmc.merged, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3)
FeaturePlot(pbmc.merged, c("CCR7"), cols.use = c("lightgrey","blue"))
FeaturePlot(pbmc.merged, c("S100A4"), cols.use = c("lightgrey","blue"))

pbmc.merged <- FindClusters(pbmc.merged, pc.use = 1:20, resolution = 1, reuse.SNN = T, do.sparse = T, print.output = F)
TSNEPlot(pbmc.merged, do.label = F)
pbmc.merged <- FindClusters(pbmc.merged, pc.use = 1:20, resolution = 0.25, reuse.SNN = T, do.sparse = T, print.output = F)
pbmc.merged <- BuildClusterTree(pbmc.merged, do.reorder = T, reorder.numeric = T)
TSNEPlot(pbmc.merged, do.label = F)



