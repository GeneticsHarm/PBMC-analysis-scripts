##
## Script for detecting outliers in a t-SNE plot
##

load("../pbmc6k/data/seurat_33k_merged")

library(Seurat)
library(pbapply)
library(FNN)


tsne.neighbours <- function(data, max.distance = 1, cores=4) {
  pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)
  for (i in 1:nrow(data)) {
    # Calculate neighbours within max.distance
    neighbours <- data$tSNE_1 > (data[i,"tSNE_1"] - max.distance) &
      data$tSNE_1 < (data[i,"tSNE_1"] + max.distance) &
      data$tSNE_2 > (data[i,"tSNE_2"] - max.distance) &
      data$tSNE_2 < (data[i,"tSNE_2"] + max.distance)
    
    data[i,"neighbours"] <- sum(neighbours) -1 # -1 to exclude current cell
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(data)
}

system.time(pbmc33k.tsne.rot <- tsne.neighbours(data = pbmc33k.merged@tsne.rot, max.distance = 1))
pbmc33k.tsne.rot[1:10,]

plot(pbmc33k.merged@tsne.rot$tSNE_2~pbmc33k.merged@tsne.rot$tSNE_1,
     pch=19, cex=0.4, col = ifelse(pbmc33k.tsne.rot$neighbours < 20 ,'red','green'))

pbmc33k.tsne.rot[1:5,]
pbmc33k.merged@snn.sparse[1:200,1:2]
## WORKSPACE

snn.sum <- rowSums(as.matrix(pbmc33k.merged@snn.sparse))

plot(pbmc33k.tsne.rot$neighbours~snn.sum, pch = 19, cex = 0.4,  col = ifelse(pbmc33k.tsne.rot$neighbours < 10 ,'red','green'))

t.cell.data <- FetchData(pbmc33k.merged, c("nUMI", "nGene", "percent.mito"), cells.use = WhichCells(pbmc33k.merged, ident = 13))

t.cell.outliers.data <- FetchData(pbmc33k.merged, c("nUMI", "nGene", "percent.mito"),
                         cells.use = WhichCells(pbmc33k.merged, cells.use = pbmc33k.tsne.rot$neighbours < 50, ident = 13))
all.cells <- FetchData(pbmc33k.merged, c("nUMI", "nGene", "percent.mito"))

outlier.data <- FetchData(pbmc33k.merged, c("nUMI", "nGene", "percent.mito"),
                                  cells.use = WhichCells(pbmc33k.merged, cells.use = pbmc33k.tsne.rot$neighbours < 10))



hist(outlier.data$nUMI, breaks = 100)
hist(all.cells$nUMI, breaks = 100)

boxplot(outlier.data$nGene, all.cells$nGene)

row.1.not.null.snn <- pbmc33k.merged@snn.sparse[1,pbmc33k.merged@snn.sparse[1,] != 0]
row.1000.not.null.snn <- pbmc33k.merged@snn.sparse[1000,pbmc33k.merged@snn.sparse[1000,] != 0]
not.null.values.snn <- pbmc33k.merged@snn.sparse


# rowSums(your.matrix != 0)

pbmc33k.merged@snn.k
dim(pbmc33k.merged@pca.rot)

pbmc33k.merged@pca.rot[1:10,1:10]

k.param <- 30
k.scale <- 25
k <- k.param*k.scale #750

genes.use <- pbmc33k.merged@var.genes
data.use <- as.matrix(pbmc33k.merged@pca.rot[, 1:25])

my.knn <- get.knn(data = data.use, k = k)
my.knn$nn.dist[1:10,1:10]
my.knn$nn.index[1:10,1:10] #Euclidian distance

dim(my.knn$nn.index)
dim(my.knn$nn.dist)

n.cells <- nrow(data.use)
nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
nn.large <- my.knn$nn.index
dim(nn.ranked)
nn.ranked[1:10,1:10]
nn.large[1:10,1:10]

length(row.1.not.null.snn)
length(row.1000.not.null.snn)

head(t.cell.data)

hist(t.cell.data$nUMI, breaks = 200)
hist(t.cell.outliers.data$nUMI, breaks = 200)

mean(t.cell.outliers.data$nUMI)
mean(t.cell.data$nUMI)
boxplot(t.cell.outliers.data$percent.mito, t.cell.data$percent.mito)

