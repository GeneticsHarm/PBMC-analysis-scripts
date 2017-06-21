####
##
## Script to merge and analyse the Pilot 3 Data
##
####

library(Seurat)
library(Matrix)
library(Matrix.utils)
library(RColorBrewer)
library(dplyr)

## Load the data

add.data <- function(orig.data = NULL, lane) {
  data <- readMM(paste0("~/Documents/Pilot 3/data/lane_", lane, "/matrix.mtx"))
  rownames(data) <- sapply(readLines(paste0("~/Documents/Pilot 3/data/lane_", lane, "/genes.tsv")), extract_field, 1, delim = "\\t")
  colnames(data) <- paste0(sapply(readLines(paste0("~/Documents/Pilot 3/data/lane_", lane, "/barcodes.tsv")), extract_field, 1, delim = "-"), "_lane", lane)
  
  if (is.null(orig.data)) return(data)
  else return(merge.Matrix(orig.data, data, by.x=rownames(orig.data), by.y=rownames(data), all.x=T, all.y=T))
}

pilot3.merged <- add.data(lane = 1)
pilot3.merged <- add.data(pilot3.merged, lane = 2)
pilot3.merged <- add.data(pilot3.merged, lane = 3)
pilot3.merged <- add.data(pilot3.merged, lane = 4)
pilot3.merged <- add.data(pilot3.merged, lane = 5)
pilot3.merged <- add.data(pilot3.merged, lane = 6)
pilot3.merged <- add.data(pilot3.merged, lane = 7)
pilot3.merged <- add.data(pilot3.merged, lane = 8)

##
## Combined clustering
##

combined.seurat <- new("seurat", raw.data = pilot3.merged)
# No cut-offs
combined.seurat <- Setup(combined.seurat, min.cells = 0, min.genes = 0, project = "pilot3", do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")

#save(combined.seurat, file = "./data/pilot3_centered_scaled.Rda")
#load(file = "./data/pilot3_centered_scaled.Rda")

## Load deAnonymized data

load("./data/pilot3_deAnonymize.Rda")
combined.seurat <- AddMetaData(combined.seurat, pilot3.de.anonymize)

## Mitochondrial genes

genes <- read.table("~/Documents/Pilot 3/data/lane_1/genes.tsv")
mito.genes <- genes[grep("^MT-", genes$V2),"V1"]

combined.percent.mito <- colSums(expm1(combined.seurat@data[mito.genes, ])) / colSums(expm1(combined.seurat@data))
combined.seurat <- AddMetaData(combined.seurat, combined.percent.mito, "percent.mito")

## QC Gene/read & mitochondrial percentage

combined.meta.data <- FetchData(combined.seurat, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
# Gene / UMI hexagon plot
ggplot(combined.meta.data, aes(nUMI, nGene)) + 
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Number of genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 2500, colour="red")  +
  geom_hline(yintercept = 500, colour="red")
# Gene / UMI doublets colored
plot(combined.meta.data$nGene~combined.meta.data$nUMI, 
     ylab = "Number of genes", xlab = "Number of reads", 
     col = ifelse(combined.meta.data$llkDoublet.llkSinglet > 100 & !combined.meta.data$lane %in% c(2,3) ,'red', rgb(0,0,0,0.1)),
     cex=0.2, cex.lab=1.4, pch=20)
# UMI / % mito genes plot
ggplot(combined.meta.data, aes(nUMI, percent.mito)) + geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Fraction mitochondrial genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 0.05, colour="red")

## Cut-offs

dim(combined.seurat@data) # 28855 cells
combined.seurat.subset <- SubsetData(combined.seurat, subset.name = "percent.mito", accept.high = 0.05)
dim(combined.seurat.subset@data) # 27709 cells
combined.seurat.subset <- SubsetData(combined.seurat.subset, subset.name = "nGene", accept.low = 500)
dim(combined.seurat.subset@data) # 28583 cells -272  (on original data) 27474 cells - 235 (on mito subsetted data)
combined.seurat.subset <- SubsetData(combined.seurat.subset, subset.name = "llkDoublet.llkSinglet", accept.high = 100)
dim(combined.seurat.subset@data) # 25342 cells (losing two samples)

sd.cut <- median(combined.meta.data$nGene) + sd(combined.meta.data$nGene) * 3
mad.cut <- median(combined.meta.data$nGene) + mad(combined.meta.data$nGene) * 3
sd.min.cut <- median(combined.meta.data$nGene) - sd(combined.meta.data$nGene) * 3
mad.min.cut <- median(combined.meta.data$nGene) - mad(combined.meta.data$nGene) * 3
sum(combined.meta.data$nGene > sd.cut) / nrow(combined.meta.data)
sum(combined.meta.data$nGene > mad.cut) / nrow(combined.meta.data)
sum(combined.meta.data$nGene < sd.min.cut) / nrow(combined.meta.data)
sum(combined.meta.data$nGene < mad.min.cut) / nrow(combined.meta.data)
sum(combined.meta.data$nGene < 500) / nrow(combined.meta.data)

combined.seurat.subset.sd <- SubsetData(combined.seurat.subset, subset.name = "nGene", accept.high = sd.cut)
dim(combined.seurat.subset.sd@data) # 24791 cells

ggplot(combined.meta.data, aes(nUMI, nGene)) + 
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Number of genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = sd.cut, colour="red")  +
  annotate("text", x=30000, y=sd.cut+100, label="> Median + 3*SD = 2.1%", size=5) +
  geom_hline(yintercept = mad.cut, colour="black") +
  annotate("text", x=30000, y=mad.cut+100, label="> Median + 3*MAD = 6,7%", size=5) +
  geom_hline(yintercept = sd.min.cut, colour="red", linetype=2)  +
  annotate("text", x=30000, y=sd.min.cut+100, label="< Median - 3*SD = 0%", size=5) +
  geom_hline(yintercept = mad.min.cut, colour="black", linetype=2) +
  annotate("text", x=30000, y=mad.min.cut+100, label="< Median - 3*MAD = 0,7%", size=5) 

ggplot(combined.meta.data, aes(nUMI, nGene)) + 
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Number of genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = sd.cut, colour="red")  +
  annotate("text", x=30000, y=sd.cut+100, label="> Median + 3*SD = 2.1%", size=5) +
  geom_hline(yintercept = 500, colour="red") +
  annotate("text", x=30000, y=500+100, label="< 500 = 0,9%", size=5)

save(combined.seurat.subset, file = "./data/pilot3_ensembl_deanonymized.Rda")

## Regress out latent variables Variable genes
combined.seurat.subset <- MeanVarPlot(combined.seurat.subset, x.low.cutoff = 0, y.cutoff = 1.2, do.plot = F)
length(combined.seurat.subset@var.genes)
combined.seurat.subset <- RegressOut(combined.seurat.subset, genes.regress = combined.seurat.subset@var.genes, latent.vars = c("percent.mito", "nUMI"))
#combined.seurat.subset <- RegressOut(combined.seurat.subset, latent.vars = c("percent.mito", "nUMI"))
#combined.seurat.subset <- RegressOut(combined.seurat.subset, genes.regress = combined.seurat.subset@var.genes, latent.vars = c("percent.mito", "nUMI", "lane"))

colnames(combined.seurat.subset@data)

counts <- colSums(combined.seurat.subset@raw.data[,colnames(combined.seurat.subset@data)])
mean(counts)

dim(combined.seurat.subset@raw.data[,colnames(combined.seurat.subset@data)])
colMeans()

## PCA/t-SNE & Clustering
combined.seurat.subset <- PCAFast(combined.seurat.subset, pc.genes = combined.seurat.subset@var.genes, pcs.compute = 40, pcs.print = 30)
PCElbowPlot(combined.seurat.subset, num.pc = 40)

combined.seurat.subset <- RunTSNE(combined.seurat.subset, dims.use = 1:16, do.fast = T)
combined.seurat.subset <- FindClusters(combined.seurat.subset, pc.use = 1:25, resolution = 1, save.SNN = T, do.sparse = T)
combined.seurat.subset <- BuildClusterTree(combined.seurat.subset, do.reorder = T, reorder.numeric = T)

#load(file = "./data/pilot3_clustered_ensembl_deanonymized.Rda")
#load(file = "./, data/pilot3_clustered_deanonymized_lane_regress.Rda")
load(file = "./data/pilot3_clustered_deanonymized.Rda")
#load(file = "./data/pilot3_clustered_deanonymized_lane_regress.Rda")
#load(file = "./data/pilot3_clustered_ensembl_deanonymized.Rda")


combined.meta.data <- FetchData(combined.seurat.subset, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
true.doublet <- combined.meta.data$llkDoublet.llkSinglet > 100 & !combined.meta.data$lane %in% c(2,3)
  
tsne.plot <- ggplot(combined.seurat.subset@tsne.rot, aes(x=tSNE_1,y=tSNE_2, colour=true.doublet)) +
  geom_point(alpha = 0.5, size=0.8) +
  theme_minimal(base_size = 14) +
  ylab("t-SNE 1") + xlab("t-SNE 2") +
  scale_color_manual(values=c("lightgrey",'red')) +
  theme(legend.title=element_blank(), legend.position="none")
  #theme(panel.background = element_rect(fill = "#f6f6f7"))
tsne.plot

doublets <- combined.meta.data[true.doublet,]
singlets <- combined.meta.data[!true.doublet,]

boxplot(doublets$nGene, singlets$nGene, outline = F, col = c("red", "lightgrey"), names = c("Doublet", "Singlet"), ylab = "nGene", cex.lab=1.2)
boxplot(doublets$nUMI, singlets$nUMI, outline = F, col = c("red", "lightgrey"), names = c("Doublet", "Singlet"), ylab = "nUMI", cex.lab=1.2)

##
## Plots
##

# function for gathering violin plot data
get.violin.data <- function(seurat, genes) {
  output <- data.frame(gene = character(0), value= numeric(0), ident = character(0))
  for (gene in genes) {
    data.use = data.frame(FetchData(seurat,gene))
    data.use = t(data.use)
    data.melt=data.frame(rep(gene, length(seurat@ident)))
    colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data.use[1,1:length(seurat@ident)])
    data.melt$id=names(data.use)[1:length(seurat@ident)]
    data.melt$ident=seurat@ident
    noise <- rnorm(length(data.melt$value))/100000
    data.melt$value=as.numeric(as.character(data.melt$value))+noise
    output <- rbind(output, data.melt)
  }
  return(output)
}

plot.data <- get.violin.data(combined.seurat.subset, c("CD3E", "CD8A", "MS4A1", "CD79A", "CD34", "NKG7", "GNLY", "CD14", "FCGR3A", "LYZ", "FCER1A", "CST3", "ITGAX", "CD1C", "PPBP"))
plot.data <- get.violin.data(combined.seurat.subset, c("FCGR3A", "NCAM1", "NCAM1-AS1", "NCAM2", "CD3D", "CD3E", "PRF1", "GZMB", "NCR1", "NCR2", "NCR3", "KLRF1", "CD226", "KLRC1", "LILRB1"))
plot.data <- get.violin.data(combined.seurat.subset, c("CD1D", "KLRB1", "LAG3", "CD27", "CD28", "B3GAT1", "ZBTB16", "CD33", "ITGAX", "CD1C", "CLEC4C", "FOXP3", "STAT1", "STAT5A", "STAT5B"))
plot.data <- get.violin.data(combined.seurat.subset, c("SELL", "LYN", "ITGAM", "ITGAX", "CX3CR1", "CSF1R", "IFITM1", "IFITM2", "IFITM3", "S100A9", "CSF3R", "GFRA2", "CLEC10A", "CD8A", "CD8B"))

violin.plot <- ggplot(plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() + 
  facet_wrap(~ gene,scales = "free_x", ncol = length(levels(plot.data$gene))) + 
  theme(strip.text.x = element_text(size=10, angle=50),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_x_discrete(limits = rev(levels(plot.data$ident))) 
# scale_fill_manual(values = colors.x.all)

tsne.plot <- ggplot(data.plot, aes(x=tSNE_1,y=tSNE_2, colour=ident)) +
  geom_point(alpha = 0.5, size=0.8) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white")) +
  # scale_color_manual(values=colors.x.all) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha=1))) + 
  theme(panel.background = element_rect(fill = "#f6f6f7"))

tsne.plot

#tsne.plot
source("./multiplot.R")
multiplot(tsne.plot, violin.plot, cols = 2)

PlotClusterTree(combined.seurat.subset)
FeaturePlot(combined.seurat.subset, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3)

## Color by lane
combined.meta.data <- FetchData(combined.seurat.subset, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane"))
ggplot(combined.seurat.subset@tsne.rot, aes(x=tSNE_1,y=tSNE_2, colour=as.factor(combined.meta.data$lane))) +
  geom_point(alpha = 1, size=0.3) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white")) +
  scale_color_manual(values=brewer.pal(8, "Set1")) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha=1))) + 
  theme(panel.background = element_rect(fill = "#f6f6f7"))

##
## Merge to major cell types
##
PlotClusterTree(combined.seurat.subset)

combined.major <- combined.seurat.subset

#levels(combined.major@ident) <- c(levels(combined.major@ident), "Erythrocyte","cMonocyte","ncMonocyte","mDC","Megakaryocyte","CD56(dim) NK", "CD56(bright) NK","pDC","Plasma","B","HPC", "CD4+ T", "CD8+ T")
levels(combined.major@ident) <- c(levels(combined.major@ident), "Monocyte","Erythrocyte","ncMonocyte","DC","Megakaryocyte","NK", "CD56(bright) NK","pDC","Plasma","B","HPC", "CD4+ T", "CD8+ T")
combined.major@ident[combined.major@ident %in% c(9,16,18:21)] <- "CD4+ T"
combined.major@ident[combined.major@ident %in% c(6,10,11,17)] <- "CD8+ T"
combined.major@ident[combined.major@ident == 1] <- "Erythrocyte"
#combined.major@ident[combined.major@ident == 2] <- "cMonocyte"
combined.major@ident[combined.major@ident %in% c(2,3)] <- "Monocyte"
#combined.major@ident[combined.major@ident == 3] <- "ncMonocyte"
#combined.major@ident[combined.major@ident == 4] <- "mDC"
combined.major@ident[combined.major@ident == 5] <- "Megakaryocyte"
#combined.major@ident[combined.major@ident == 7] <- "CD56(dim) NK"
combined.major@ident[combined.major@ident %in% c(7,8)] <- "NK"
#combined.major@ident[combined.major@ident == 8] <- "CD56(bright) NK"
#combined.major@ident[combined.major@ident == 12] <- "pDC"
combined.major@ident[combined.major@ident  %in% c(12, 4)] <- "DC"
combined.major@ident[combined.major@ident == 13] <- "Plasma"
combined.major@ident[combined.major@ident == 14] <- "B"
combined.major@ident[combined.major@ident == 15] <- "HPC"





TSNEPlot(combined.major)

colors.temp = c("#153057", "#009ddb","#965ec8","#e64b50","#edba1b", "orange", "#71bc4b", "black", "darkgrey", "lightgrey", "pink", "red", "white")
#colors.x = c("#e64b50", "#71bc4b", "#965ec8", "#edba1b", "#009ddb", "black", "lightgrey", "#153057")
#colors.temp <- colorRampPalette(brewer.pal(9, "Set1"))(13)


get.percentage <- function(seurat, ident) {
  paste0(round(sum(seurat@ident == ident) / length(seurat@ident),3) * 100, "% ", ident)
}

create.labels <- function(seurat, idents) {
  sapply(rev(levels(idents)), get.percentage, seurat = seurat)
}

violin.plot.data <- get.violin.data(combined.major, c("CD3D", "CD8A","NKG7","GZMB","KLRC1","CD14","LYZ","FCGR3A","LYN","MS4A1","CD79A","ITGAX","CD1C","CLEC4C","PPBP","CD34"))
#violin.plot.data <- get.violin.data(combined.major, c("CD3E", "CD8A", "MS4A1", "CD79A", "NKG7", "GNLY", "CD14", "FCGR3A", "LYZ", "FCER1A", "CST3", "ITGAX", "CD34", "PPBP"))

violin.plot.subset <- ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() + 
  facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) + 
  theme(strip.text.x = element_text(size=12, angle=50),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,0,0,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top", labels=create.labels(combined.major, violin.plot.data$ident)) +
  scale_fill_manual(values = colors.temp)
violin.plot.subset

violin.plot.data <- get.violin.data(combined.major, c("CD3D", "CD8A","GZMB","KLRC1","CD14","FCGR3A","MS4A1","CD79A","ITGAX","PPBP","CD34"))

violin.plot.subset <- ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F, lwd=0) + 
  ylab("") + xlab("") +
  coord_flip() + 
  theme_minimal() +
  facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) + 
  theme(strip.text.x = element_text(size=14, angle=50),
        axis.text.y = element_text(size=10),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin = unit(c(0,0,0,1.3), "cm"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank()) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), labels=NULL) +
  scale_fill_manual(values = colors.temp)
violin.plot.subset

data.plot <- combined.major@tsne.rot
data.plot$ident <- combined.major@ident
data.plot %>% dplyr::group_by(ident) %>% summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) -> centers
centers$tSNE_1 <- centers$tSNE_1 - 3

tsne.plot.subset <- ggplot(data.plot, aes(x=tSNE_1,y=tSNE_2, colour=ident)) +
  theme_minimal() +
  geom_point(alpha = 0.5, size=0.8) +
  #theme(panel.background = element_rect(fill = "white")) +
  scale_color_manual(values=colors.temp, labels=paste0(1:13,": ", levels(data.plot$ident))) +
  theme(legend.title=element_blank(), legend.text=element_text(size=14)) +
  guides(colour = guide_legend(override.aes = list(size=8, alpha=1))) + 
  #theme(panel.background = element_rect(fill = "#f6f6f7")) + 
  ylab("t-SNE 1") + xlab("t-SNE 2") +
  geom_point(data = centers, aes(x=tSNE_1, y=tSNE_2), size=8, shape=21, fill=rgb(1,1,1,0.7), color="black") + 
  geom_text(data=centers, aes(label=1:length(levels(ident))), size = 4, color="black")
tsne.plot.subset


layout <- matrix(c(1,1,1,1,2,2), ncol = 2, nrow = 3, byrow=TRUE)
multiplot(tsne.plot.subset, violin.plot.subset, layout = layout)



## Markers

unknown.markers <- FindMarkers(combined.major, "Unknown")
dendritic.markers <- FindMarkers(combined.major, "Dendritic cells")
all.markers <- FindAllMarkers(combined.seurat.subset, max.cells.per.ident = 100)
write.table(unknown.markers, file = "./data/markers_unkown.tsv", sep = "\t", quote = F)
write.table(dendritic.markers, file = "./data/dendritic.markers.tsv", sep = "\t", quote = F)
write.table(all.markers, file = "./data/all.markers.tsv", sep = "\t", quote = F)
grep(pattern = "CD63", rownames(combined.seurat.subset@data), value = T)

TSNEPlot(combined.seurat.subset, do.label = T)

## Color by sample


combined.major.meta.data <- FetchData(combined.major, c("nUMI", "nGene", "percent.mito", "llkDoublet.llkSinglet", "lane", "assigned_sample.s."))

ggplot(combined.major@tsne.rot, aes(x=tSNE_1,y=tSNE_2, colour=combined.major.meta.data$assigned_sample.s.)) +
  geom_point(alpha = 0.5, size=0.8) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha=1))) + 
  theme(panel.background = element_rect(fill = "#f6f6f7")) + 
  ylab("t-SNE 1") + xlab("t-SNE 2")

plot(combined.major@tsne.rot$tSNE_1, combined.major@tsne.rot$tSNE_2, 
     ylab = "t-SNE 1", xlab = "t-SNE 2", 
     col = ifelse(combined.major.meta.data$assigned_sample.s. == "1_LLDeep_0198",'red', rgb(0,0,0,0.1)),
     cex=0.2, cex.lab=1.4, pch=20)
plot(combined.major@tsne.rot$tSNE_1, combined.major@tsne.rot$tSNE_2, 
     ylab = "t-SNE 1", xlab = "t-SNE 2", 
     col = ifelse(combined.major.meta.data$lane == 8,'red', rgb(0,0,0,0.1)),
     cex=0.2, cex.lab=1.4, pch=20)


get.percentage.temp <- function(ident, sample) {
   sum(sample@ident == ident) / length(sample@ident) * 100
}

get.percentage.temp <- function(ident, sample) {
  sum(sample@ident == ident)
}

create.percentages <- function(seurat, sample.ids) {
  samples.matrix <- matrix(0L, nrow = length(levels(sample.ids)), ncol = length(levels(seurat@ident)), dimnames = list(levels(sample.ids), levels(seurat@ident)))
  
  for (sample in levels(sample.ids)) {
    sample.subset <- SubsetData(seurat, cells.use = WhichCells(seurat, cells.use = sample.ids == sample))
    percentages <- sapply(levels(seurat@ident), get.percentage.temp, sample = sample.subset)
    samples.matrix[sample,] <- percentages
  }
  return(samples.matrix)
}

sample.percentages <- create.percentages(seurat = combined.major, sample.ids = combined.major.meta.data$assigned_sample.s.)

mv <- c("M","V","M","V","M","V","M","M","M","M","V","M","M","V","V","M","M","V","M","V","M","V","V","M","M","V","V","V","V","M","V","V","M","M","M","M","V","M","M","V","V","M","M","V","V")

sum(mv == "V")
sum(mv == "M")

sample.percentages.male <- sample.percentages[mv == "M",]
sample.percentages.female <- sample.percentages[mv == "V",]

boxplot(sample.percentages.male[,1], sample.percentages.female[,1],
        sample.percentages.male[,2], sample.percentages.female[,2],
        sample.percentages.male[,3], sample.percentages.female[,3],
        sample.percentages.male[,4], sample.percentages.female[,4],
        sample.percentages.male[,5], sample.percentages.female[,5],
        sample.percentages.male[,6], sample.percentages.female[,6],
        sample.percentages.male[,7], sample.percentages.female[,7],
        sample.percentages.male[,8], sample.percentages.female[,8],
        sample.percentages.male[,9], sample.percentages.female[,9],
        sample.percentages.male[,10], sample.percentages.female[,10],
        sample.percentages.male[,11], sample.percentages.female[,11],
        sample.percentages.male[,12], sample.percentages.female[,12],
        names = c("Th male", "Th female",
                  "Tc male", "Tc female",
                  "T? male", "T? female",
                  "NKT? male", "NKT? female",
                  "NK male", "NK female",
                  "B male", "B female",
                  "Plasma M", "Plasma V",
                  "Mono CD14+ M", "Mono CD14+ V",
                  "Mono FCGR3A+ M", "Mono FCGR3A+ V",
                  "Dendritic M", "Dendritic V",
                  "HSCs M", "HSCs V",
                  "Megakaryo M", "Megakaryo V"
        ), las =2)


par(mar=c(10,5,1,1))

sample.percentages

require(reshape2)
melt(sample.percentages)

colors = c("#153057","#009ddb","#e64b50","#edba1b","#71bc4b","#965ec8")

ggplot(data = melt(sample.percentages), aes(x=Var2, y=value)) + 
  geom_boxplot(aes(fill=Var2), lwd=0.5) +
  theme_minimal(base_family = "Helvetica Neue", base_size = 16) +
  scale_fill_manual(breaks=colors, values = colors) +
  coord_cartesian(ylim = c(0, 500)) +
  ylab("Cells per donor") +
  xlab("")

ggplot(sample.percentages, aes(x=snp, y=expression, group=snp)) +
  geom_boxplot(notch=F, color = "black", outlier.shape=NA, fill= colors, lwd=0.6, alpha=1) + 
  theme_minimal(base_family = "Helvetica Neue Light") +
  theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica Neue")) +
  #theme_bw(base_family = 'Helvetica') +
  #theme(panel.background = element_rect(fill = "#f6f6f7")) +
  facet_wrap(~cell.type, ncol = length(levels(data$cell.type)) ) +
  geom_point(position = position_jitter(width = 0.2), size = 0.8, color="black", alpha = 0.5) +
  ggtitle(paste0(gene, " (", snp.name, ")")) +
  ylab("Expression") +
  xlab("")

boxplot(sample.percentages, las=2)
boxplot(sample.percentages[grep("gonl",rownames(sample.percentages)),], las=2)
boxplot(sample.percentages[-grep("gonl",rownames(sample.percentages)),], las=2)

write.table(sample.percentages, file = "./data/pilot3_percentages_samples.tsv", sep = "\t", quote = F)

## Plot al the lanes 
plot.list <- list()
for (i in 1:8) {
  local({
    i <- i
    plot = ggplot(combined.major@tsne.rot) + 
    geom_point(aes(x = tSNE_1, y = tSNE_2, colour = combined.major.meta.data$lane == i), size=0.2) +
    scale_colour_manual(values = c(rgb(0,0,0,0.1), 'red')) +
    theme(legend.position="none")
    plot.list[[i]] <<- plot
  })
}
multiplot(plotlist = plot.list, cols = 4)

##
## Create eQTL expression file 
##

create.samples <- function(exp.data, sample.ids) {
  #exp.data <- exp.data[apply(exp.data, 1, function(x) !all(x==0)),]
  
  samples.matrix <- matrix(0L, nrow = nrow(exp.data), ncol = length(levels(sample.ids)), dimnames = list(row.names(exp.data), levels(sample.ids)))

  samples.totals <- 0
  
  for (sample in levels(sample.ids)) {
    print(sample)
    sample.cells <- sample.ids == sample
    
    samples.totals <- samples.totals + sum(sample.cells)
    
    print(sum(sample.cells))
    if (sum(sample.cells) <= 1) {
      samples.matrix[,sample] <- exp.data[,sample.cells]
    } else {
      samples.matrix[,sample] <- rowMeans(exp.data[,sample.cells]) # Remember to load MatrixUtils
    }
  }
  
  print(samples.totals / ncol(samples.matrix))
  
  return(samples.matrix)
}


## All cells
sample.ids <- FetchData(combined.seurat.subset, "assigned_sample.s.")$assigned_sample.s.
#all.cells.exp.data <- combined.seurat.subset@data[apply(combined.seurat.subset@data, 1, function(x) !all(x==0)),]
all.cells.exp.data <- combined.seurat.subset@data
dim(all.cells.exp.data)
all.cells.samples <- create.samples(all.cells.exp.data, sample.ids)
all.cells.samples <- create.samples(combined.seurat.subset@scale.data, sample.ids)
write.table(all.cells.samples, file = "./data/pilot3_traitfile_ensembl_all_log_norm.tsv", sep = "\t", quote = F)
write.table(all.cells.samples, file = "/groups/umcg-wijmenga/tmp02/projects/scRNAseq_10X_pilot/eQTL/pilot3_traitfile_ensembl_all_scaled.tsv", sep = "\t", quote = F)
head(all.cells.samples)

##
## Get ensemble ids 
##


## Subset the Th cells
combined.th.cells <- SubsetData(combined.major, ident.use = "CD4+ T")
combined.th.cells <- SubsetData(combined.seurat.subset, ident.use = "CD4+ T cells")
combined.th.cells <- RegressOut(combined.th.cells, latent.vars = c("percent.mito", "nUMI"))

sample.ids <- FetchData(combined.th.cells, "assigned_sample.s.")$assigned_sample.s.

combined.th.cells@scale.data[1:10,1:10]
## Create the sample / mean gene expression matrix & save it
dim(combined.th.cells@scale.data)
#th.exp.data <- combined.th.cells@scale.data[apply(combined.th.cells@scale.data, 1, function(x) !all(x==0)),]
th.exp.data <- combined.th.cells@scale.data
dim(th.exp.data)
th.cells.samples <- create.samples(th.exp.data, sample.ids)
write.table(th.cells.samples, file = "./data/pilot3_traitfile_th_scaled_centered.tsv", sep = "\t", quote = F)


th.cells <- SubsetData(combined.seurat.subset, ident.use = c(9,16,18:21))
sample.ids <- FetchData(th.cells, "assigned_sample.s.")$assigned_sample.s.
th.cells.samples <- create.samples(th.cells@data, sample.ids)
write.table(th.cells.samples, file = "./data/pilot3_traitfile_th_log_norm.tsv", sep = "\t", quote = F)

b.cells <- SubsetData(combined.seurat.subset, ident.use = 14)

## B cells
b.cells <- SubsetData(combined.seurat.subset, ident.use = 14)
b.cells <- RegressOut(b.cells, latent.vars = c("percent.mito", "nUMI"))
#b.exp.data <- b.cells@scale.data[apply(b.cells@scale.data, 1, function(x) !all(x==0)),]
b.exp.data <- b.cells@scale.data
sample.ids <- FetchData(b.cells, "assigned_sample.s.")$assigned_sample.s.
b.cells.samples <- create.samples(b.exp.data, sample.ids)
dim(b.cells.samples)
write.table(b.cells.samples, file = "./data/pilot3_traitfile_b_scaled.tsv", sep = "\t", quote = F)

b.cells.samples <- create.samples(b.cells@data, sample.ids)
write.table(b.cells.samples, file = "./data/pilot3_traitfile_b_log_norm.tsv", sep = "\t", quote = F)


## NK cells
nk.cells <- SubsetData(combined.seurat.subset, ident.use = c(7,8))
nk.cells <- RegressOut(nk.cells, latent.vars = c("percent.mito", "nUMI"))
#nk.exp.data <- nk.cells@scale.data[apply(nk.cells@scale.data, 1, function(x) !all(x==0)),]
nk.exp.data <- nk.cells@scale.data
sample.ids <- FetchData(nk.cells, "assigned_sample.s.")$assigned_sample.s.
nk.cells.samples <- create.samples(nk.exp.data, sample.ids)
dim(nk.cells.samples)
write.table(nk.cells.samples, file = "./data/pilot3_traitfile_nk_scaled.tsv", sep = "\t", quote = F)

nk.cells.samples <- create.samples(nk.cells@data, sample.ids)
write.table(nk.cells.samples, file = "./data/pilot3_traitfile_nk_log_norm.tsv", sep = "\t", quote = F)


## Monocytes 
monocytes <- SubsetData(combined.seurat.subset, ident.use = c(2,3))
monocytes <- RegressOut(monocytes, latent.vars = c("percent.mito", "nUMI"))
monocytes.exp.data <- monocytes@scale.data[apply(monocytes@scale.data, 1, function(x) !all(x==0)),]
monocytes.exp.data <- monocytes@scale.data
sample.ids <- FetchData(monocytes, "assigned_sample.s.")$assigned_sample.s.
monocytes.samples <- create.samples(monocytes.exp.data, sample.ids)
dim(monocytes.samples)
write.table(monocytes.samples, file = "./data/pilot3_traitfile_mono_scaled.tsv", sep = "\t", quote = F)

monocytes.samples <- create.samples(monocytes@data, sample.ids)
write.table(monocytes.samples, file = "./data/pilot3_traitfile_mono_log_norm.tsv", sep = "\t", quote = F)


## CD8+ 
tc.cells <- SubsetData(combined.seurat.subset, ident.use = c(6,10,11,17))
tc.cells <- RegressOut(tc.cells, latent.vars = c("percent.mito", "nUMI"))
tc.exp.data <- tc.cells@scale.data[apply(tc.cells@scale.data, 1, function(x) !all(x==0)),]
tc.exp.data <- tc.cells@scale.data
sample.ids <- FetchData(tc.cells, "assigned_sample.s.")$assigned_sample.s.
tc.samples <- create.samples(tc.exp.data, sample.ids)
dim(tc.samples)
write.table(tc.samples, file = "./data/pilot3_traitfile_tc_scaled.tsv", sep = "\t", quote = F)

tc.samples <- create.samples(tc.cells@data, sample.ids)
write.table(tc.samples, file = "./data/pilot3_traitfile_tc_log_norm.tsv", sep = "\t", quote = F)

## Dendritic cells 
dend.cells <- SubsetData(combined.seurat.subset, ident.use = c(4,12))
dend.cells <- RegressOut(dend.cells, latent.vars = c("percent.mito", "nUMI"))
dend.exp.data <- dend.cells@scale.data[apply(dend.cells@scale.data, 1, function(x) !all(x==0)),]
dend.exp.data <- dend.cells@scale.data
sample.ids <- FetchData(dend.cells, "assigned_sample.s.")$assigned_sample.s.
dend.samples <- create.samples(dend.exp.data, sample.ids)
dim(dend.samples)
write.table(dend.samples, file = "./data/pilot3_traitfile_dend_scaled.tsv", sep = "\t", quote = F)

dend.samples <- create.samples(dend.cells@data, sample.ids)
write.table(dend.samples, file = "./data/pilot3_traitfile_dend_log_norm.tsv", sep = "\t", quote = F)


VlnPlot(combined.major, "CD3D")

## Check data quality
th.cells.samples.sorted <- t(th.cells.samples)
th.cells.samples.sorted <- th.cells.samples.sorted[,order(colSums(th.cells.samples.sorted), decreasing = T )] # Order on totals

boxplot(th.cells.samples.sorted[,1:10], outline=T)
boxplot(th.cells.samples.sorted[,500:510], outline=T)
boxplot(th.cells.samples.sorted[,1000:1010], outline=T)
boxplot(th.cells.samples.sorted[,2000:2010], outline=T)
boxplot(th.cells.samples.sorted[,3000:3010], outline=T)
boxplot(th.cells.samples.sorted[,10000:10010], outline=T)

##
## Try a few different parameters for clustering
##

DoClusterAnalysis <- function(seurat, pc.use=25, y.cutoff = 0.8, ssn.res = 0.25) {
  seurat <- MeanVarPlot(seurat, x.low.cutoff = 0, y.cutoff = y.cutoff, do.plot = F)
  seurat <- PCAFast(seurat, pc.genes = seurat@var.genes, pcs.compute = pc.use, do.print = F)
  PCElbowPlot(seurat, num.pc = pc.use)
  seurat <- RunTSNE(seurat, dims.use = 1:pc.use, do.fast = T)
  #seurat <- FindClusters(seurat, pc.use = 1:pc.use, resolution = ssn.res, save.SNN = T, do.sparse = T)
  #seurat <- BuildClusterTree(seurat, do.reorder = T, reorder.numeric = T)
  FeaturePlot(seurat, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3)
  return(seurat)
}

length(combined.seurat@var.genes)
length(combined.seurat.y8.pc25@var.genes)
length(combined.seurat.y10.pc25@var.genes)

combined.seurat.y8.pc16 <- DoClusterAnalysis(combined.seurat, pc.use = 16, y.cutoff = 0.8)
combined.seurat.y8.pc20 <- DoClusterAnalysis(combined.seurat, pc.use = 20, y.cutoff = 0.8)
combined.seurat.y8.pc25 <- DoClusterAnalysis(combined.seurat, pc.use = 25, y.cutoff = 0.8)

combined.seurat.y10.pc16 <- DoClusterAnalysis(combined.seurat, pc.use = 16, y.cutoff = 1)
combined.seurat.y10.pc20 <- DoClusterAnalysis(combined.seurat, pc.use = 20, y.cutoff = 1)
combined.seurat.y10.pc25 <- DoClusterAnalysis(combined.seurat, pc.use = 25, y.cutoff = 1)

combined.seurat.y12.pc14 <- DoClusterAnalysis(combined.seurat, pc.use = 14, y.cutoff = 1.2)
combined.seurat.y12.pc16 <- DoClusterAnalysis(combined.seurat, pc.use = 16, y.cutoff = 1.2)
combined.seurat.y12.pc20 <- DoClusterAnalysis(combined.seurat, pc.use = 20, y.cutoff = 1.2)
combined.seurat.y12.pc25 <- DoClusterAnalysis(combined.seurat, pc.use = 25, y.cutoff = 1.2)


combined.seurat.y14.pc16 <- DoClusterAnalysis(combined.seurat, pc.use = 16, y.cutoff = 1.4)


FeaturePlot(combined.seurat.y8.pc20, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y8.pc16, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y8.pc25, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y10.pc25, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y10.pc16, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y12.pc16, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y12.pc14, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y12.pc20, c("GNLY"), cols.use = c("lightgrey","blue"))
FeaturePlot(combined.seurat.y14.pc16, c("GNLY"), cols.use = c("lightgrey","blue"))

plot.feature <- function(seurat, feature, title) {
  palette <- colorRampPalette(c("lightgrey", "blue","orange","red"))
  nUMI <- FetchData(seurat, feature)
  colors <- palette(20)[as.numeric(cut(nUMI,breaks = 20))]
  plot(seurat@tsne.rot$tSNE_1, seurat@tsne.rot$tSNE_2, cex = 0.4, pch=16, col = colors,
       main = title, ylab = "", xlab = "")
}

par(mfrow=c(3,3), mar=c(1,1,3,1))
plot.feature(combined.seurat.y8.pc25, "GNLY", "Variable cut-off=0.8 PC=25")
plot.feature(combined.seurat.y8.pc20, "GNLY", "Variable cut-off=0.8 PC=20")
plot.feature(combined.seurat.y8.pc16, "GNLY", "Variable cut-off=0.8 PC=16")
plot.feature(combined.seurat.y10.pc25, "GNLY", "Variable cut-off=1.0 PC=25")
plot.feature(combined.seurat.y10.pc20, "GNLY", "Variable cut-off=1.0 PC=20")
plot.feature(combined.seurat.y10.pc16, "GNLY", "Variable cut-off=1.0 PC=16")
plot.feature(combined.seurat.y12.pc25, "GNLY", "Variable cut-off=1.2 PC=25")
plot.feature(combined.seurat.y12.pc20, "GNLY", "Variable cut-off=1.2 PC=20")
plot.feature(combined.seurat.y12.pc16, "GNLY", "Variable cut-off=1.2 PC=16")
par(mfrow=c(1,1))

par(mfrow=c(3,3), mar=c(1,1,3,1))
plot.feature(combined.seurat.y8.pc25, "NKG7", "Variable cut-off=0.8 PC=25")
plot.feature(combined.seurat.y8.pc20, "NKG7", "Variable cut-off=0.8 PC=20")
plot.feature(combined.seurat.y8.pc16, "NKG7", "Variable cut-off=0.8 PC=16")
plot.feature(combined.seurat.y10.pc25, "NKG7", "Variable cut-off=1.0 PC=25")
plot.feature(combined.seurat.y10.pc20, "NKG7", "Variable cut-off=1.0 PC=20")
plot.feature(combined.seurat.y10.pc16, "NKG7", "Variable cut-off=1.0 PC=16")
plot.feature(combined.seurat.y12.pc25, "NKG7", "Variable cut-off=1.2 PC=25")
plot.feature(combined.seurat.y12.pc20, "NKG7", "Variable cut-off=1.2 PC=20")
plot.feature(combined.seurat.y12.pc16, "NKG7", "Variable cut-off=1.2 PC=16")
par(mfrow=c(1,1))

par(mfrow=c(3,3), mar=c(1,1,3,1))
plot.feature(combined.seurat.y8.pc25, "CD8A", "Variable cut-off=0.8 PC=25")
plot.feature(combined.seurat.y8.pc20, "CD8A", "Variable cut-off=0.8 PC=20")
plot.feature(combined.seurat.y8.pc16, "CD8A", "Variable cut-off=0.8 PC=16")
plot.feature(combined.seurat.y10.pc25, "CD8A", "Variable cut-off=1.0 PC=25")
plot.feature(combined.seurat.y10.pc20, "CD8A", "Variable cut-off=1.0 PC=20")
plot.feature(combined.seurat.y10.pc16, "CD8A", "Variable cut-off=1.0 PC=16")
plot.feature(combined.seurat.y12.pc25, "CD8A", "Variable cut-off=1.2 PC=25")
plot.feature(combined.seurat.y12.pc20, "CD8A", "Variable cut-off=1.2 PC=20")
plot.feature(combined.seurat.y12.pc16, "CD8A", "Variable cut-off=1.2 PC=16")
par(mfrow=c(1,1))


FeaturePlot(combined.seurat.y10.pc16, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3)
FeaturePlot(combined.seurat.y12.pc16, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3)
FeaturePlot(combined.seurat.y14.pc16, c("MS4A1", "GNLY","CD3E","CD8A","LYZ","PF4"), cols.use = c("lightgrey","blue"), nCol = 3)

## Sparsity plot

rowSums(combined.seurat.subset@data)

non.zero.counts <- rowSums(as.matrix(combined.seurat.subset@data) != 0) # per gene aantal non zero counts
non.zero.frac <- non.zero.counts / ncol(combined.seurat.subset@data)

hist(non.zero.frac, breaks = 100, ylim = c(1, 5000),
     col = "red", ylab = "Frequency of genes", xlab = "Fraction of non-zero values",
     main = "", cex.lab = 1.4)

combined.seurat.subset@data[1:3,1:3]


