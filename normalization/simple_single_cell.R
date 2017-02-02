##
## Test of the SSC (simpleSingleCell) normalization method
## This method uses a pooling approach
##

#source("http://bioconductor.org/workflows.R")
#workflowInstall("simpleSingleCell")

#install.packages("limSolve")
library(limSolve)
library(scater)
library(scran)
library(Seurat)

#load("data/simple_single_cell.RData")
load("data/seurat_33k_merged")

#rename cluster IDs
new_ids=c("Megakaryocytes","Monocytes",3,"NK_1","NK_2","B_Cells_Plasma","Dendritic_Cells","B_Cells_1","B_Cells_2","B_Cells_3","T_Cells_CD4+_1","T_Cells_CD8+_1","T_Cells_CD4+_2","T_Cells_CD8+_2","T_Cells_CD8+_3","T_Cells_CD8+_4")
current.cluster.ids=1:16
pbmc33k.merged@ident <- plyr::mapvalues(pbmc33k.merged@ident, from = current.cluster.ids, to = new_ids)

nk.cells <- pbmc33k.merged@raw.data[,WhichCells(pbmc33k.merged, c("NK_1", "NK_2"))]
nk.cells <- nk.cells[apply(nk.cells[,-1], 1, function(x) !all(x==0)),]

b.cells <- pbmc33k.merged@raw.data[,WhichCells(pbmc33k.merged, c("B_Cells_1", "B_Cells_2", "B_Cells_3"))]
b.cells <- b.cells[apply(b.cells[,-1], 1, function(x) !all(x==0)),]

t.h.cells <- pbmc33k.merged@raw.data[,WhichCells(pbmc33k.merged, c("T_Cells_CD4+_1", "T_Cells_CD4+_2"))]
t.h.cells <- t.h.cells[apply(t.h.cells[,-1], 1, function(x) !all(x==0)),]
dim(t.h.cells)

sce.t.h.cells <- newSCESet(countData = as.matrix(t.h.cells))
dim(sce.t.h.cells)

is.mito <- grepl("^MT\\.", rownames(sce.t.h.cells))
sce.t.h.cells <- calculateQCMetrics(sce.t.h.cells, feature_controls=list(Mt=is.mito))


par(mfrow=c(1,2))
hist(sce.t.h.cells$total_counts/1e3, xlab="Library sizes (thousands)", main="", breaks=20, col="grey80", ylab="Number of cells")
hist(sce.t.h.cells$total_features, xlab="Number of expressed genes", main="", breaks=20, col="grey80", ylab="Number of cells")
par(mfrow=c(1,1))

#plot(sizeFactors(sce.t.h.cells), sce.t.h.cells$total_counts/1e3, log="xy", ylab="Library size (thousands)", xlab="Size factor")

# We remove cells with log-library sizes that are more than 3 median absolute deviations (MADs)
# below the median log-library size. (A log-transformation improves resolution at small values, especially when the MAD of the raw values is comparable to or greater than the median.)
# We also remove cells where the log-transformed number of expressed genes is 3 MADs below the median.
libsize.drop <- isOutlier(sce.t.h.cells$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce.t.h.cells$total_features, nmads=3, type="lower", log=TRUE)

hist(sce.t.h.cells$pct_counts_feature_controls_Mt, xlab="Mitochondrial proportion (%)",  ylab="Number of cells", breaks=20, main="", col="grey80")

dim(sce.t.h.cells)
sce.t.h.cells <- sce.t.h.cells[,!(libsize.drop | feature.drop)]
dim(sce.t.h.cells)

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce.t.h.cells, pca_data_input="pdata") + fontsize

#anno <- select(org.Mm.eg.db, keys=rownames(sce.t.h.cells), keytype="SYMBOL", column="ENSEMBL")
#ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]
#assignments <- cyclone(sce, mm.pairs, gene.names=ensembl)
#plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

ave.counts <- rowMeans(counts(sce.t.h.cells))
keep <- rowMeans(counts(sce.t.h.cells)) >= 0.05
sum(keep)

hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=log10(0.05), col="blue", lwd=2, lty=2)

#plotQC(sce, type = "highest-expression", n=50)

numcells <- nexprs(sce.t.h.cells, byrow=TRUE)
smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), ylab="Number of expressing cells")

sce.t.h.cells <- sce.t.h.cells[keep,]

sce.t.h.cells <- sce.t.h.cells[!fData(sce.t.h.cells)$is_feature_control_Mt,]

# Any systematic difference in count size across the non-DE majority of genes between two cells is assumed to represent bias 
# and is removed by scaling. More specifically, “size factors” are calculated that represent the extent to which counts should be scaled
# in each library.
#b.cells.clusters <- quickCluster(sce.t.h.cells)
sce.t.h.cells <- computeSumFactors(sce.t.h.cells, positive=TRUE)
#sce <- computeSumFactors(sce, sizes=c(20, 40, 60, 80), positive=TRUE)
summary(sizeFactors(sce.t.h.cells))
sce.t.h.cells <- sce.t.h.cells[,sizeFactors(sce.t.h.cells) != 0]
#dim(sce.t.h.cells)
plot(sizeFactors(sce.t.h.cells), sce.t.h.cells$total_counts/1e3, log="xy", ylab="Library size (thousands)", xlab="Size factor")
# Lun ATL, Bach K and Marioni JC (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biol. 17:75

# ?? misschien niet doen
sce.t.h.cells <- normalize(sce.t.h.cells)

#plotExplanatoryVariables(sce, variables=c("pct_counts_feature_controls_Mt", "total_counts"))

var.fit <- trendVar(sce.t.h.cells, trend="loess", use.spikes=FALSE, span=0.4)
var.out <- decomposeVar(sce.t.h.cells, var.fit)

plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
#points(var.fit$mean, var.fit$var, col="red", pch=16)
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="red", lwd=2)

hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.02),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
nrow(hvg.out)

head(hvg.out)


plotExpression(sce.t.h.cells, rownames(hvg.out)[1:10])

set.seed(100)
var.cor <- correlatePairs(sce.t.h.cells, subset.row=rownames(hvg.out))
#var.cor <- correlatePairs(sce, subset.row=names(expressed_genes[expressed_genes]))
head(var.cor)

var.cor[1996990:1997001,]

sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)







