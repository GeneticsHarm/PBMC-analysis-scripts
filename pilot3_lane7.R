
library(Seurat)
library(Matrix)

lane7.data <- Read10X("~/Documents/Pilot 3/data/lane_7/")
lane7.seurat  <- new("seurat", raw.data = lane7.data)
lane8.data <- Read10X("~/Documents/Pilot 3/data/lane_8/")
lane8.seurat  <- new("seurat", raw.data = lane8.data)

pbmc8k.data <- Read10X("~/Documents/Stage/Data/pbmc8k/outs/filtered_gene_bc_matrices/GRCh38/")
pbmc8k.reads.per.cell <- colSums(pbmc8k.data)
mean(pbmc8k.reads.per.cell)

lane7.data[1:40,1:40]

total.reads.per.cell <- colSums(lane8.data)
mean(total.reads.per.cell)
median(total.reads.per.cell)

sum(total.reads.per.cell)

hist(total.reads.per.cell, breaks = 40)
dim(lane7.data)

lane7.raw <- readMM("~/Documents/Pilot 3/data/lane_7/matrix.mtx")
rownames(lane7.raw) <- read.delim("~/Documents/Pilot 3/data/lane_7/genes.tsv", sep = "\n", header = F)[,1]
colnames(lane7.raw) <- read.delim("~/Documents/Pilot 3/data/lane_7/barcodes.tsv", sep = "\n", header = F)[,1]
lane7.raw.cell.totals <- colSums(lane7.raw)
mean(lane7.raw.cell.totals)
median(lane7.raw.cell.totals)

lane7.seurat  <- new("seurat", raw.data = lane7.data)
lane7.seurat <- Setup(lane7.seurat, min.cells = 3, min.genes = 200, project = "10X_PBMC", do.scale = F, do.center = F, names.field = 2, names.delim = "\\-")

lane7.mito.genes <- grep("^MT-", rownames(lane7.seurat@data), value = T)
lane7.percent.mito <- colSums(expm1(lane7.seurat@data[lane7.mito.genes, ])) / colSums(expm1(lane7.seurat@data))

lane7.seurat <- AddMetaData(lane7.seurat, lane7.percent.mito, "percent.mito")

GenePlot(lane7.seurat, "nUMI", "nGene")
all.cells <- FetchData(lane7.seurat, c("nUMI", "nGene", "percent.mito"))
plot(all.cells$nGene~all.cells$nUMI, ylab = "Number of genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.3), cex=0.2, cex.lab=1.4)
abline(h = 2500, lwd = 2, lty = 2, col = "red")
abline(h = 500, lwd = 2, lty = 2, col = "red")

plot(all.cells$percent.mito~all.cells$nUMI, ylab = "Percentage mitochondrial genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.4), cex=0.5, cex.lab=1.4)
abline(h = 0.05, lwd = 2, lty = 2, col = "red")

## !!! WERKT NIET
lane7.seurat <- SubsetData(lane7.seurat, subset.name = "nGene", accept.high = 2500)
lane7.seurat <- SubsetData(lane7.seurat, subset.name = "percent.mito", accept.high = 0.05)
lane7.seurat <- SubsetData(lane7.seurat, subset.name = "nGene", accept.low = 500)

lane7.seurat <- MeanVarPlot(lane7.seurat, x.low.cutoff = 0, y.cutoff = 0.8)
length(lane7.seurat@var.genes)

lane7.seurat <- RegressOut(lane7.seurat, latent.vars = c("percent.mito", "nUMI"), genes.regress = lane7.seurat@var.genes)
lane7.seurat <- PCAFast(lane7.seurat, pc.genes = lane7.seurat@var.genes, pcs.compute = 40, pcs.print = 30)
PCElbowPlot(lane7.seurat, num.pc = 40)

lane7.seurat <- RunTSNE(lane7.seurat, dims.use = 1:25, do.fast = T)

lane7.seurat <- FindClusters(lane7.seurat, pc.use = 1:25, resolution = seq(0.25), save.SNN = T, do.sparse = T)

TSNEPlot(lane7.seurat, do.label = F)

