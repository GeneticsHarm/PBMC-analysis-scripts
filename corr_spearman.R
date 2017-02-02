##
## Calculate the Spearman's rank correlation coefficient for every PBMC cell type
##
## Dataset: 33K PBMC 10xGenomics
## Loaded data: data/seurat_33k_merged (created @ pbmc33k.R)
## Saved data: data/33k_correlation_spearman_celltype
##

library(Matrix)
library(scran)
library(Seurat)

## Data loading

# Load the clustered data.
load("data/seurat_33k_merged")

# Rename cluster IDs.
new_ids=c("Megakaryocytes","Monocytes",3,"NK_1","NK_2","B_Cells_Plasma","Dendritic_Cells","B_Cells_1","B_Cells_2","B_Cells_3","T_Cells_CD4+_1","T_Cells_CD8+_1","T_Cells_CD4+_2","T_Cells_CD8+_2","T_Cells_CD8+_3","T_Cells_CD8+_4")
current.cluster.ids=1:16
pbmc33k.merged@ident <- plyr::mapvalues(pbmc33k.merged@ident, from = current.cluster.ids, to = new_ids)

## Data preperation

# Filter genes with less than 10 percent non-zero values.
fraction.non.zero <- 0.05
pbmc33k.non.zero.counts <- rowSums(as.matrix(pbmc33k.merged@data) != 0)
genes.filtered <- pbmc33k.non.zero.counts > fraction.non.zero * ncol(pbmc33k.merged@data)
pbmc33k.filtered.counts <- pbmc33k.merged@data[genes.filtered,]

dim(pbmc33k.merged@data)
dim(pbmc33k.filtered.counts)

#TSNEPlot(pbmc33k.merged)

#fraction.cell.type = 1/2
 
# Subset the celltypes from filtered data
nk.cells <- pbmc33k.filtered.counts[,WhichCells(pbmc33k.merged, c("NK_1", "NK_2"))]
b.cells <- pbmc33k.filtered.counts[,WhichCells(pbmc33k.merged, c("B_Cells_1", "B_Cells_2", "B_Cells_3", "B_Cells_Plasma"))]
t.h.cells <- pbmc33k.filtered.counts[,WhichCells(pbmc33k.merged, c("T_Cells_CD4+_1", "T_Cells_CD4+_2"))]
t.c.cells <- pbmc33k.filtered.counts[,WhichCells(pbmc33k.merged, c("T_Cells_CD8+_1", "T_Cells_CD8+_2", "T_Cells_CD8+_3","T_Cells_CD8+_4"))]
monocytes <- pbmc33k.filtered.counts[,WhichCells(pbmc33k.merged, c("Monocytes"))]

## Calculate correlations

# When all sample are used algorithm fails
sample.count <- 1500

nk.cells.cor.pairs <- correlatePairs(as.matrix(nk.cells)[,1:sample.count])
nk.cells.cor.sign <- nk.cells.cor.pairs[nk.cells.cor.pairs$FDR < 0.05,]
b.cells.cor.pairs <- correlatePairs(as.matrix(b.cells)[,1:sample.count])
b.cells.cor.pairs.sign <- b.cells.cor.pairs[b.cells.cor.pairs$FDR < 0.05,]
t.h.cells.cor.pairs <- correlatePairs(as.matrix(t.h.cells)[,1:sample.count])
t.h.cells.cor.pairs.sign <- t.h.cells.cor.pairs[t.h.cells.cor.pairs$FDR < 0.05,]
t.c.cells.cor.pairs <- correlatePairs(as.matrix(t.c.cells)[,1:sample.count])
t.c.cells.cor.pairs.sign <- t.c.cells.cor.pairs[t.c.cells.cor.pairs$FDR < 0.05,]
monocytes.cor.pairs <- correlatePairs(as.matrix(monocytes)[,1:sample.count])
monocytes.cor.pairs.sign <- monocytes.cor.pairs[monocytes.cor.pairs$FDR < 0.05,]

## Save the correlations

save(nk.cells.cor.sign, nk.cells.cor.pairs, b.cells.cor.pairs.sign, b.cells.cor.pairs, t.c.cells.cor.pairs.sign, t.c.cells.cor.pairs, t.h.cells.cor.pairs.sign, t.h.cells.cor.pairs, monocytes.cor.pairs.sign, monocytes.cor.pairs, file = "data/33k_correlation_spearman_celltype_frac_0.05")



