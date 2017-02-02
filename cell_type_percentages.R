##
## Script for calculation of the fraction of every cell type in the 33K dataset
##
## Loaded data: ./data/pbmc33k.Rda
##

load("./data/pbmc33k.Rda")

TSNEPlot(pbmc33k)

b.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("BCell_A", "BCell_B", "BCell_C"))]
t.h.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("CD4_Naive", "CD4_Memory_Antiviral", "CD4_Memory_TIGIT", "CD4_Memory"))]
t.c.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("CD8Effector_GZMH", "CD8Effector_GZMK", "CD8_Memory", "CD8_Naive"))]
nk.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("NK", "NK_Prss57", "NK_XCL1"))]
monocytes <- pbmc33k@data[,WhichCells(pbmc33k, c("Mono_CD16+", "Mono_CD16+_C1qa", "Mono_Apobec3b", "Mono_Apobec3a", "Mono_CD14+", "Mono_CD14+_Inflam", "Mono_CD14+_Antiviral"))]
megakaryos <- pbmc33k@data[,WhichCells(pbmc33k, "Megakaryocyte")]
dendritic.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("DC_CD141+", "DC_CD1C+"))]

total.count <- length(pbmc33k@cell.names) 

b.cell.percentage <- dim(b.cells)[2] / total.count * 100
t.h.cell.percentage <- dim(t.h.cells)[2] / total.count * 100
t.c.cell.percentage <- dim(t.c.cells)[2] / total.count * 100
nk.cells.percentage <- dim(nk.cells)[2] / total.count * 100
monocytes.percentage <- dim(monocytes)[2] / total.count * 100
megakaryos.percentage <- dim(megakaryos)[2] / total.count * 100
dendritic.cells.percentage <- dim(dendritic.cells)[2] / total.count * 100

b.cell.percentage
t.h.cell.percentage
t.c.cell.percentage
nk.cells.percentage
monocytes.percentage
megakaryos.percentage
dendritic.cells.percentage
