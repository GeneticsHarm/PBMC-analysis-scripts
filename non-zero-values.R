##
## Scipt for exploration of the zero values.
##
## Loaded data: data/seurat_33k_merged
##

# Load the clustered data.
load("data/seurat_33k_merged")

# Rename cluster IDs.
new_ids=c("Megakaryocytes","Monocytes",3,"NK_1","NK_2","B_Cells_Plasma","Dendritic_Cells","B_Cells_1","B_Cells_2","B_Cells_3","T_Cells_CD4+_1","T_Cells_CD8+_1","T_Cells_CD4+_2","T_Cells_CD8+_2","T_Cells_CD8+_3","T_Cells_CD8+_4")
current.cluster.ids=1:16
pbmc33k.merged@ident <- plyr::mapvalues(pbmc33k.merged@ident, from = current.cluster.ids, to = new_ids)

pbmc33k.merged@data[1:10,1:10]
dim(pbmc33k.merged@data)
ncol(pbmc33k.merged)

pbmc33k.non.zero.counts <- rowSums(as.matrix(pbmc33k.merged@data) != 0) # per gene aantal non zero counts
pbmc33k.non.zero.frac <- pbmc33k.non.zero.counts / ncol(pbmc33k.merged@data) # per gene fractie non zero counts

genes.filtered.5 <- pbmc33k.non.zero.frac > 0.05
genes.frac.5 <- pbmc33k.non.zero.frac[genes.filtered.5]
pbmc33k.filtered.n.zero.5 <- pbmc33k.merged@data[genes.filtered.5,]

genes.filtered.10 <- pbmc33k.non.zero.counts > 0.10 * ncol(pbmc33k.merged@data)
genes.frac.10 <- pbmc33k.non.zero.frac[genes.filtered.10]
pbmc33k.filtered.n.zero.10 <- pbmc33k.merged@data[genes.filtered.10,]

genes.filtered.15 <- pbmc33k.non.zero.counts > 0.15 * ncol(pbmc33k.merged@data)
genes.frac.15 <- pbmc33k.non.zero.frac[genes.filtered.15]
pbmc33k.filtered.n.zero.15 <- pbmc33k.merged@data[genes.filtered.15,]

dim(pbmc33k.filtered.n.zero.5)
dim(pbmc33k.filtered.n.zero.10)
dim(pbmc33k.filtered.n.zero.15)

hist(pbmc33k.non.zero.frac, breaks = 300, main = "", col=rgb(0.242, 0.746, 0.591, 1), ylab = "Number of genes", xlab = "Fraction non-zero values", cex.lab=1.4, xlim = c(0,0.4))

hist(genes.frac.5, xlim = c(0,1), col = "green", breaks = 100, main = "", ylab = "Number of genes", xlab = "Fraction non-zero values", cex.lab=1.4)
abline(v = 0.05, lwd = 2, lty = 2, col = "red")
abline(v = 0.10, lwd = 2, lty = 2, col = "red")
abline(v = 0.15, lwd = 2, lty = 2, col = "red")
box()



hist(pbmc33k.non.zero.frac, breaks = 200, xlim = c(0,0.2))

h = hist(pbmc33k.non.zero.frac)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)

plot(density(pbmc33k.non.zero.frac))
plot(density(genes.frac.10))

genes.filtered <- pbmc33k.non.zero.counts > fraction.non.zero * ncol(pbmc33k.merged@data)
pbmc33k.filtered.counts <- pbmc33k.merged@data[genes.filtered,]

factors <- factor( c(rep(0.05, times = length(genes.frac.5)), rep(0.10, times = length(genes.frac.10)), rep(0.15, times = length(genes.frac.15)) )  )

genes.frac <- data.frame(fraction = factors, non.zero = c(genes.frac.5, genes.frac.10, genes.frac.15) )
ggplot(genes.frac, aes(x=non.zero, fill=fraction)) + geom_histogram(binwidth=.01, alpha=.6, position="identity")
ggplot(genes.frac, aes(x=non.zero, fill=fraction)) + geom_density(alpha=.5)
ggplot(genes.frac, aes(x=non.zero, fill=fraction)) + geom_histogram(binwidth=.03, position="dodge")

## ggplot

library(ggplot2)

set.seed(1234)
dat <- data.frame(cond = factor(rep(c("A","B"), each=200)), rating = c(rnorm(200),rnorm(200, mean=.8)))

head(dat)

# Overlaid histograms
ggplot(dat, aes(x=rating, fill=cond)) + geom_histogram(binwidth=.5, alpha=.5, position="identity")



