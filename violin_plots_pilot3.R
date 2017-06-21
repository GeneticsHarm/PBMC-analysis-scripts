library(Seurat)

source("./multiplot.R")
load(file = "./data/pilot3_clustered_deanonymized.Rda")

markers.six <- FindMarkers(combined.seurat.subset, ident.1 = 6, ident.2 = c(9,10,11), max.cells.per.ident = 100)
markers.nine <- FindMarkers(combined.seurat.subset, ident.1 = 9, ident.2 = c(6,10,11), max.cells.per.ident = 100)
markers.ten <- FindMarkers(combined.seurat.subset, ident.1 = 10, ident.2 = c(6,9,11), max.cells.per.ident = 100)
markers.eleven <- FindMarkers(combined.seurat.subset, ident.1 = 11, ident.2 = c(6,9,10), max.cells.per.ident = 100)
write.table(markers.six, file = "./data/markers-6__9-10-11.tsv", sep = "\t", quote = F)
write.table(markers.nine, file = "./data/markers-9__6-10-11.tsv", sep = "\t", quote = F)
write.table(markers.ten, file = "./data/markers-10__6-9-11.tsv", sep = "\t", quote = F)
write.table(markers.eleven, file = "./data/markers-11__6-9-10.tsv", sep = "\t", quote = F)


markers.17 <- FindMarkers(combined.seurat.subset, ident.1 = 17, ident.2 = c(18,19,20,21), max.cells.per.ident = 100)
markers.18 <- FindMarkers(combined.seurat.subset, ident.1 = 18, ident.2 = c(17,19,20,21), max.cells.per.ident = 100)
markers.19 <- FindMarkers(combined.seurat.subset, ident.1 = 19, ident.2 = c(17,18,20,21), max.cells.per.ident = 100)
markers.20 <- FindMarkers(combined.seurat.subset, ident.1 = 20, ident.2 = c(17,18,19,21), max.cells.per.ident = 100)
markers.21 <- FindMarkers(combined.seurat.subset, ident.1 = 21, ident.2 = c(17,18,19,20), max.cells.per.ident = 100)

markers.16 <- FindMarkers(combined.seurat.subset, ident.1 = 16, ident.2 = c(17,18,19,20,21), max.cells.per.ident = 100)

write.table(markers.16, file = "./data/markers-16__17-18-19-20-21.tsv", sep = "\t", quote = F)
write.table(markers.17, file = "./data/markers-17__18-19-20-21.tsv", sep = "\t", quote = F)
write.table(markers.18, file = "./data/markers-18__17-19-20-21.tsv", sep = "\t", quote = F)
write.table(markers.19, file = "./data/markers-19__17-18-20-21.tsv", sep = "\t", quote = F)
write.table(markers.20, file = "./data/markers-20__17-18-19-21.tsv", sep = "\t", quote = F)
write.table(markers.21, file = "./data/markers-21__17-18-19-20.tsv", sep = "\t", quote = F)


th.nk.subset <- SubsetData(combined.seurat.subset, ident.use = c(6,9:11,16:21))

plot.data <- get.violin.data(th.nk.subset, c("SELL", "CCR7", "CXCR4", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4"))
grep("CD3", rownames(th.nk.subset@data), value = T)

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

tsne.plot <- ggplot(th.nk.subset@tsne.rot, aes(x=tSNE_1,y=tSNE_2, colour=th.nk.subset@ident)) +
  geom_point(alpha = 0.5, size=0.8) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white")) +
  # scale_color_manual(values=colors.x.all) +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=6, alpha=1))) + 
  theme(panel.background = element_rect(fill = "#f6f6f7"))

multiplot(tsne.plot, violin.plot, cols = 2)

plot.data <- get.violin.data(th.nk.subset, c("STAT1", "STAT4", "EOMES", "TBX21", "PRDM1", "PDCD1", "HAVCR2", "SOX9", "LAG3", "FOXO1", "FOXO3", "IL10", "TGFB1", "IL12A", "EBI3", "IL2RA", "CTLA4", "TNFRSF18", "IL4", "IFNG", "CCR7", "RUNX3", "IFNGR1", "CXCR3", "IL4R", "IL1RL1", "CCR4", "IL17RB", "PTGDR2", "GATA3", "STAT6", "BHLHE41", "MAF", "SPI1", "IL23R", "CCR6",  "IL1R1", "RORA", "RORC", "CCR10", "TNFSF4", "CD40LG", "ICOS", 'PDCD1', "BCL6", "STAT3", "BCL6B", "MBD2", "BMI1"))

grep("IL1R", rownames(th.nk.subset@data), value = T)
