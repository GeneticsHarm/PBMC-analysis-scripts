load("data/seurat_33k_main")

TSNEPlot(pbmc33k.merged, do.label = T)

pbmc33k.subset <- SubsetData(pbmc33k.merged, ident.use = c("B", "NK", "TH", "TC", "Monocytes"))
TSNEPlot(pbmc33k.subset, do.label = T)

markers <- FindAllMarkers(pbmc33k.subset, do.print = T)
FeaturePlot(pbmc33k.subset, "ENSG00000161570")

markers[markers$cluster=="TH",][1:10,]
top.markers <- NULL
for (cluster in levels(pbmc33k.subset@ident)) {
    # select the first 10 markers where the percentage in the other cluster of that marker is < 10%
    cluster.markers <- markers[markers$cluster==cluster & markers$pct.2<0.10,]
    top.markers <- rbind(top.markers, cluster.markers[1:10,])
}
top.markers

doublet.scores <- matrix(nrow = length(pbmc33k.subset@cell.names),
                         ncol = length(levels(pbmc33k.subset@ident)) + 1,
                         dimnames = list(pbmc33k.subset@cell.names, c("highest.frac", levels(pbmc33k.subset@ident))))


idents <- levels(pbmc33k.subset@ident)
pb <- txtProgressBar(min = 0, max = nrow(pbmc33k.subset@data), style = 3)
for (i in 1:nrow(pbmc33k.subset@data)) {
  
  cell <- pbmc33k.subset@data[,i]
  for (ident in idents) {
    top.markers.ident <- top.markers[top.markers$cluster == ident, "gene"]
    cell.top.markers.ident <- cell[top.markers.ident]
    doublet.scores[i, ident] <- sum(cell.top.markers.ident)
  }
  doublet.scores[i,1] <- max(doublet.scores[i,], na.rm = T) / sum(doublet.scores[i,], na.rm = T) 
  setTxtProgressBar(pb, i)
}
close(pb)

hist(doublet.scores[,1])
doublet.scores[1:25,]

head(doublet.scores[is.na(doublet.scores[,1]),])

doublet.scores[is.nan(doublet.scores)] <- 1
doublet.scores[is.na(doublet.scores)] <- 1

hist(doublet.scores[,1])
length(doublet.scores[,1])

pbmc33k.subset <- AddMetaData(pbmc33k.subset, doublet.scores[,1], "doublet.score")

FeaturePlot(pbmc33k.subset, "percent.mito", cols.use = c("blue", "lightgrey"))
all.cells <- FetchData(pbmc33k.subset, c("nUMI", "nGene", "percent.mito", "doublet.score"))
plot(all.cells$nGene, all.cells$doublet.score, pch=20, col=rgb(0,0,0,alpha=0.3), cex=0.2, cex.lab=1.4)






