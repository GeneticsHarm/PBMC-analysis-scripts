library(Seurat)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(hexbin)

##
## Load the data
##
pilot3.aggr.data <- Read10X("~/Documents/Pilot 3/data/aggr/")
pilot3.aggr.seurat <- new("seurat", raw.data = pilot3.aggr.data)
pilot3.aggr.seurat <- Setup(pilot3.aggr.seurat, min.cells = 3, min.genes = 200, project = "pilot3.aggr", do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")

##
## QC
##
pilot3.mito.genes <- grep("^MT-", rownames(pilot3.aggr.seurat@data), value = T)
pilot3.percent.mito <- colSums(expm1(pilot3.aggr.seurat@data[pilot3.mito.genes, ])) / colSums(expm1(pilot3.aggr.seurat@data))
pilot3.aggr.seurat <- AddMetaData(pilot3.aggr.seurat, pilot3.percent.mito, "percent.mito")

pilot3.aggr.meta <- FetchData(pilot3.aggr.seurat, c("nUMI", "nGene", "percent.mito", "orig.ident"))

## Gene/reads per cell plot
ggplot(pilot3.aggr.meta, aes(nUMI, nGene)) +
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Number of genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 2500, colour="red")  + 
  geom_hline(yintercept = 500, colour="red")

## Percent mito / nUMI plot
ggplot(pilot3.aggr.meta, aes(nUMI, percent.mito)) + 
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Fraction mitochondrial genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 0.05, colour="red")

# Cut-off
pilot3.aggr.seurat <- SubsetData(pilot3.aggr.seurat, subset.name = "nGene", accept.high = 2500)
pilot3.aggr.seurat <- SubsetData(pilot3.aggr.seurat, subset.name = "percent.mito", accept.high = 0.05)
pilot3.aggr.seurat <- SubsetData(pilot3.aggr.seurat, subset.name = "nGene", accept.low = 500)

##
## CLustering
##

## Variable genes 
pilot3.aggr.seurat <- MeanVarPlot(pilot3.aggr.seurat, x.low.cutoff = 0, y.cutoff = 1.2, do.plot = F)
length(pilot3.aggr.seurat@var.genes)
## Linear regression: nUMI and percentage UMI's on all the cells
pilot3.aggr.seurat <- RegressOut(pilot3.aggr.seurat, genes.regress = pilot3.aggr.seurat@var.genes, latent.vars = c("percent.mito", "nUMI"))

## PCA/t-SNE & Clustering
pilot3.aggr.seurat <- PCAFast(pilot3.aggr.seurat, pc.genes = pilot3.aggr.seurat@var.genes, pcs.compute = 40, do.print = F)
PCElbowPlot(pilot3.aggr.seurat, num.pc = 40)
pilot3.aggr.seurat <- RunTSNE(pilot3.aggr.seurat, dims.use = 1:16, do.fast = T)
pilot3.aggr.seurat <- FindClusters(pilot3.aggr.seurat, pc.use = 1:25, resolution = seq(0.25), save.SNN = T, do.sparse = T)
## Plots
pilot3.aggr.seurat <- BuildClusterTree(pilot3.aggr.seurat, do.reorder = T, reorder.numeric = T)
TSNEPlot(pilot3.aggr.seurat, do.label = F, pt.size = 0.5)
FeaturePlot(pilot3.aggr.seurat, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3, pt.size = 0.5)


##
## Save & load
##

#save(pilot3.aggr.seurat, file = "./data/pilot_3_aggr_clustered.Rda")
#load("./data/pilot_3_aggr_clustered.Rda")


##
## Lane bias
##
library(RColorBrewer)
palette(brewer.pal(8, "Set1"))

lanes <- factor(sub(".+\\-(\\d)", "Lane \\1", colnames(pilot3.aggr.seurat@data)), ordered = T)

plot(pilot3.aggr.seurat@tsne.rot$tSNE_1, pilot3.aggr.seurat@tsne.rot$tSNE_2, cex = 0.2, col = lanes, ylab = "t-SNE 2", xlab = "t-SNE 1")
legend(-35,-20,unique(lanes),col=1:length(lanes),pch=16, text.width = 5)

pilot3.aggr.per.lane <- SetIdent(pilot3.aggr.seurat, ident.use = lanes)
TSNEPlot(pilot3.aggr.per.lane, do.label = F, pt.size = 0.2, colors.use = palette())
