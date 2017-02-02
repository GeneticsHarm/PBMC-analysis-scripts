##
## Script for calution of Pearson correlation on the Seurat test data
##
## Loaded data: ./data/pbmc33k.Rda
## Saved data: ./data/cell_type_cor.Rda
##

load("./data/pbmc33k.Rda")

b.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("BCell_A", "BCell_B", "BCell_C"))]
t.h.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("CD4_Naive", "CD4_Memory_Antiviral", "CD4_Memory_TIGIT", "CD4_Memory"))]
t.c.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("CD8Effector_GZMH", "CD8Effector_GZMK", "CD8_Memory", "CD8_Naive"))]
nk.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("NK", "NK_Prss57", "NK_XCL1"))]
monocytes <- pbmc33k@data[,WhichCells(pbmc33k, c("Mono_CD16+", "Mono_CD16+_C1qa", "Mono_Apobec3b", "Mono_Apobec3a", "Mono_CD14+", "Mono_CD14+_Inflam", "Mono_CD14+_Antiviral"))]
megakaryos <- pbmc33k@data[,WhichCells(pbmc33k, "Megakaryocyte")]
dendritic.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("DC_CD141+", "DC_CD1C+"))]

dim(b.cells)
dim(t.h.cells)
dim(t.c.cells)
dim(nk.cells)
dim(monocytes)
dim(megakaryos)
dim(dendritic.cells)

b.cells <- b.cells[apply(b.cells[,-1], 1, function(x) !all(x==0)),] # remove genes with no expression
b.cells.cor <- cor(t(as.matrix(b.cells)))
nk.cells <- nk.cells[apply(nk.cells[,-1], 1, function(x) !all(x==0)),]
nk.cells.cor <- cor(t(as.matrix(nk.cells)))
monocytes <- monocytes[apply(monocytes[,-1], 1, function(x) !all(x==0)),]
monocytes.cor <- cor(t(as.matrix(monocytes)))
megakaryos <- megakaryos[apply(megakaryos[,-1], 1, function(x) !all(x==0)),]
megakaryos.cor <- cor(t(as.matrix(megakaryos)))
dendritic.cells <- dendritic.cells[apply(dendritic.cells[,-1], 1, function(x) !all(x==0)),]
dendritic.cells.cor <- cor(t(as.matrix(dendritic.cells)))
t.h.cells <- t.h.cells[apply(t.h.cells[,-1], 1, function(x) !all(x==0)),]
t.h.cells.cor <- cor(t(as.matrix(t.h.cells)))
t.c.cells <- t.c.cells[apply(t.c.cells[,-1], 1, function(x) !all(x==0)),]
t.c.cells.cor <- cor(t(as.matrix(t.c.cells)))

save(b.cells.cor, t.h.cells.cor, t.c.cells.cor, monocytes.cor, megakaryos.cor, dendritic.cells.cor, nk.cells.cor, file = "./data/cell_type_cor.Rda")


