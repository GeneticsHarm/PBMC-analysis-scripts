####
##
## Pilot 4 Clustering
##
####

library(Seurat)
library(Matrix)
library(Matrix.utils)

lane1.data <- Read10X("~/Documents/Pilot 4/data/lane_1/")
colnames(lane1.data) <- paste0(colnames(lane1.data), "_lane1")
lane2.data <- Read10X("~/Documents/Pilot 4/data/lane_2/")
colnames(lane2.data) <- paste0(colnames(lane2.data), "_lane2")
lane3.data <- Read10X("~/Documents/Pilot 4/data/lane_3/")
colnames(lane3.data) <- paste0(colnames(lane3.data), "_lane3")

pilot4.merged <- merge.Matrix(lane1.data, lane2.data, by.x=rownames(lane1.data), by.y=rownames(lane2.data), all.x=T, all.y=T)
pilot4.merged <- merge.Matrix(pilot3.merged, lane3.data, by.x=rownames(pilot3.merged), by.y=rownames(lane3.data), all.x=T, all.y=T)

lane2.seurat <- new("seurat", raw.data = lane2.data)
lane2.seurat <- Setup(lane2.seurat, min.cells = 3, min.genes = 200, project = "pilot4.lane2", do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")

lane2.mito.genes <- grep("^MT-", rownames(lane2.seurat@data), value = T)
lane2.percent.mito <- colSums(expm1(lane2.seurat@data[lane2.mito.genes, ])) / colSums(expm1(lane2.seurat@data))
lane2.seurat <- AddMetaData(lane2.seurat, lane2.percent.mito, "percent.mito")

## Gene/read & mitochondrial plot

all.cells <- FetchData(lane2.seurat, c("nUMI", "nGene", "percent.mito"))
plot(all.cells$nGene~all.cells$nUMI, ylab = "Number of genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.3), cex=0.2, cex.lab=1.4)
abline(h = 2500, lwd = 2, lty = 2, col = "red")
abline(h = 500, lwd = 2, lty = 2, col = "red")

plot(all.cells$percent.mito~all.cells$nUMI, ylab = "Percentage mitochondrial genes", xlab = "Number of reads", pch=20, col=rgb(0,0,0,alpha=0.4), cex=0.5, cex.lab=1.4)
abline(h = 0.05, lwd = 2, lty = 2, col = "red")

## Cut-offs

lane2.seurat <- SubsetData(lane2.seurat, subset.name = "nGene", accept.high = 2500)
lane2.seurat <- SubsetData(lane2.seurat, subset.name = "percent.mito", accept.high = 0.05)
lane2.seurat <- SubsetData(lane2.seurat, subset.name = "nGene", accept.low = 500)
lane2.seurat <- RegressOut(lane2.seurat, latent.vars = c("percent.mito", "nUMI"))

## Variable genes

lane2.seurat <- MeanVarPlot(lane2.seurat, x.low.cutoff = 0, y.cutoff = 0.6, do.plot = F)

## PCA/t-SNE & Clustering

lane2.seurat <- PCAFast(lane2.seurat, pc.genes = lane2.seurat@var.genes, pcs.compute = 40, pcs.print = 30)
PCElbowPlot(lane2.seurat, num.pc = 40)
lane2.seurat <- RunTSNE(lane2.seurat, dims.use = 1:17, do.fast = T)
lane2.seurat <- lane2.seurat(lane2.seurat, pc.use = 1:20, resolution = seq(0.25), save.SNN = T, do.sparse = T)
lane2.seurat <- BuildClusterTree(lane2.seurat, do.reorder = T, reorder.numeric = T)

## Plots

TSNEPlot(lane2.seurat, do.label = F)
FeaturePlot(lane2.seurat, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"),cols.use = c("lightgrey","blue"),nCol = 3)

## Set sample identities from deAnonymize

#lane2.xy <- read.table("~/Documents/Pilot 3/xy/lane_2.txt", header = T)
lane2.ids <- read.table("~/Documents/Pilot 4/deAnonymizer/p4|l2.txt", header = T)

rownames(lane2.ids) <- sub("([ATCG]+)-1", "\\1", rownames(lane2.ids))
rownames(lane2.ids) <- paste0(rownames(lane2.ids), "_lane2")
lane2.ids.seurat <- lane2.ids[rownames(lane2.ids) %in% colnames(lane2.seurat@data),]
lane2.seurat <- SetIdent(lane2.seurat, ident.use = lane2.ids.seurat$sample)

plot(lane2.seurat@tsne.rot$tSNE_1, lane2.seurat@tsne.rot$tSNE_2, cex = 0.5, pch=16, col = ifelse(lane2.ids.seurat$doublet == 1 ,'red','lightgrey'))

lane2.seurat@cell.names[lane2.ids.seurat$doublet == 1]

lane2.doublets <- SubsetData(lane2.seurat, cells.use = lane2.seurat@cell.names[lane2.ids.seurat$doublet == 1])

dim(lane2.doublets@data)

lane2.doublets.meta <- FetchData(lane2.doublets, c("nUMI", "nGene", "percent.mito"))
lane2.meta <- FetchData(lane2.seurat, c("nUMI", "nGene", "percent.mito"))

boxplot(lane2.doublets.meta$nGene, lane2.meta$nGene, names = c("doublet nGene", "nGene total") )
boxplot(lane2.doublets.meta$nUMI, lane2.meta$nUMI, names = c("doublet nUMI", "nUMI total") )


lane2.xy <- read.table("~/Documents/Pilot 4/xy/lane_2.txt", header = T)
lane1.xy <- read.table("~/Documents/Pilot 4/xy/lane_1.txt", header = T)
lane1.ids <- read.table("~/Documents/Pilot 4/deAnonymizer/p4|l1.txt", header = T)

lane2.xy[,"cell"] <- sub("([ATCG]+)-\\d{1}", "\\1", lane2.xy[,"cell"]) # trim the barcode
lane2.xy <- lane2.xy[order(lane2.xy$cell),] # order

lane1.xy[,"cell"] <- sub("([ATCG]+)-\\d{1}", "\\1", lane1.xy[,"cell"]) # trim the barcode
lane1.xy <- lane1.xy[order(lane1.xy$cell),] # order

dim(lane2.xy)
dim(lane2.ids)

lane1.ids$doublet == 1
boxplot(lane2.xy$y[!lane2.ids$doublet == 1]~lane2.ids$sample[!lane2.ids$doublet == 1], ylab = "Y reads", xlab = "Sample")
boxplot(lane1.xy$y[!lane1.ids$doublet == 0]~lane1.ids$sample[!lane1.ids$doublet == 0], ylab = "Y reads", xlab = "Sample")

boxplot(lane1.xy$y[!lane1.ids$doublet == 0]~lane1.ids$sample[!lane1.ids$doublet == 0], ylab = "Y reads", xlab = "Sample")
lane1.xy[1:10,]
lane1.ids[1:10,]

lane1.ids.doublet.sep <- lane1.ids
lane1.ids.doublet.sep[lane1.ids$doublet==1,]$sample <- 6
boxplot(lane1.xy$y~lane1.ids.doublet.sep$sample, ylab = "Y reads", xlab = "Sample")

lane1.combined <- cbind(lane1.xy,lane1.ids.doublet.sep)
lane1.combined[1:10,]

sum(lane1.combined[lane1.combined$doublet==1,"y"]==0) / sum(lane1.combined$doublet==1)
sum(lane1.combined[lane1.combined$doublet==1,"y"]>0)

rownames(lane2.ids) <- sub("([ATCG]+)-\\d{1}", "\\1", rownames(lane2.ids)) # trim the barcode
sum(rownames(lane2.ids) == lane2.xy[,"cell"])

lane1.exp.data <- Read10X("~/Documents/Pilot 4/data/lane_1/")

lane1.seurat <- new("seurat", raw.data = lane1.exp.data)
lane1.seurat <- Setup(lane1.seurat, min.cells = 3, min.genes = 200, project = "pilot4.lane1", do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")

lane1.mito.genes <- grep("^MT-", rownames(lane1.seurat@data), value = T)
lane1.percent.mito <- colSums(expm1(lane1.seurat@data[lane1.mito.genes, ])) / colSums(expm1(lane1.seurat@data))
lane1.seurat <- AddMetaData(lane1.seurat, lane1.percent.mito, "percent.mito")

lane1.all.cells <- FetchData(lane1.seurat, c("nUMI", "nGene", "percent.mito"))

rownames(lane1.ids)[1:10]
colnames(lane1.seurat@data)[1:10]


rownames(lane1.ids) <- sub("([ATCG]+)-1", "\\1", rownames(lane1.ids))

lane1.ids.seurat <- lane1.ids[rownames(lane1.ids) %in% colnames(lane1.seurat@data),]
lane1.seurat <- SetIdent(lane1.seurat, ident.use = lane1.ids.seurat$sample)

lane1.doublets <- SubsetData(lane1.seurat, cells.use = lane1.seurat@cell.names[lane1.ids.seurat$doublet == 1])

dim(lane1.doublets@data)

lane1.doublets.meta <- FetchData(lane1.doublets, c("nUMI", "nGene", "percent.mito"))
lane1.meta <- FetchData(lane1.seurat, c("nUMI", "nGene", "percent.mito"))

boxplot(lane1.doublets.meta$nGene, lane1.meta$nGene, names = c("doublet nGene", "nGene total"), outline = F)
boxplot(lane1.doublets.meta$nUMI, lane1.meta$nUMI, names = c("Doublets lane1 nUMI", "nUMI lane1 all cells"), outline = F)


library(RColorBrewer)
library(ggplot2)
lane1.all <- cbind(lane1.ids, lane1.xy)

give.n <- function(x) {
  return(data.frame(y = unname(quantile(x, probs = 0.75, type = 3))-0.25, label = paste0("N=",length(x))))
}

ggplot(aes(x=sample, y=y, colour=sample, group=sample), data=lane1.all) +
  geom_boxplot(notch=T, colour = brewer.pal(6, "Set1"), outlier.alpha = 0.3) + 
  stat_summary(geom = 'point', fun.y=mean, shape = 23, size = 5, colour=brewer.pal(6, "Set1")) +  
  theme_bw(base_family = 'Helvetica') +
  theme(legend.position= 'none') +
  coord_cartesian(ylim = c(0, 10)) +
  stat_summary(fun.data = give.n, geom = "text", color = "black") 

ggplot(aes(x=sample, y=y, colour=sample, group=sample), data=lane1.all) +
  geom_boxplot(notch=T, colour = colors.x, outlier.alpha = 0.3) + 
  stat_summary(geom = 'point', fun.y=mean, shape = 23, size = 5, colour=brewer.pal(6, "Set1")) +  
  theme_bw(base_family = 'Helvetica') +
  theme(legend.position= 'none') +
  coord_cartesian(ylim = c(0, 10)) +
  stat_summary(fun.data = give.n, geom = "text", color = "black") 

lane1.data[1:10,1:3]
lane1.xy[1:10,]  
cell.totals.lane1 <- colSums(lane1.data)  
cell.totals.lane1[1:10]                                                       
