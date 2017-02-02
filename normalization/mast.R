##
## Test of the MAST normalization method for scRNA-seq data
##

library(MAST)

suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
})

#options(mc.cores = detectCores() - 1) #if you have multiple cores to spin
options(mc.cores = 4)

vignette('MAST-intro')
vignette('MAITAnalysis')

load("data/seurat_33k_merged")

new_ids=c("Megakaryocytes","Monocytes",3,"NK_1","NK_2","B_Cells_Plasma","Dendritic_Cells","B_Cells_1","B_Cells_2","B_Cells_3","T_Cells_CD4+_1","T_Cells_CD8+_1","T_Cells_CD4+_2","T_Cells_CD8+_2","T_Cells_CD8+_3","T_Cells_CD8+_4")
current.cluster.ids=1:16
pbmc33k.merged@ident <- plyr::mapvalues(pbmc33k.merged@ident, from = current.cluster.ids, to = new_ids)

b.cells <- pbmc33k.merged@data[,WhichCells(pbmc33k.merged, c("B_Cells_1","B_Cells_2","B_Cells_3"))]
b.cells.mast <- FromMatrix(as.matrix(b.cells))
b.cells.mast

aheatmap(assay(b.cells.mast[1:1000,1:100]), labRow='', annCol=as.data.frame(colData(b.cells.mast)[,c('condition', 'ourfilter')]), distfun='spearman')

cdr <- colSums(assay(b.cells.mast)>0)
colData(b.cells.mast)$cngeneson <- scale(cdr)

b.cells.sample <- b.cells.mast[sample(which(freq(b.cells.mast)>.1), 20),]
b.cells.flat <- as(b.cells.sample, 'data.table')

ggplot(b.cells.flat, aes(x=value))+geom_density() + facet_wrap(~primerid, scale='free_y')

thres <- thresholdSCRNACountMatrix(assay(b.cells.mast), nbins = 20, min_per_bin = 30)
par(mfrow=c(4,4))
plot(thres)

assays(b.cells.mast) <- list(thresh=thres$counts_threshold, tpm=assay(b.cells.mast))

expressed_genes <- freq(b.cells.mast) > 0.2
b.cells.mast.expressed <- b.cells.mast[expressed_genes,]
dim(b.cells.mast.expressed)

assay(b.cells.mast.expressed[1:20,1:3])

#### MAIT

data(maits, package='MAST')
dim(maits$expressionmat)
head(maits$cdat)
dim(maits$cdat)
head(maits$fdat)

scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)
aheatmap(assay(scaRaw[1:1000,]), labRow='', annCol=as.data.frame(colData(scaRaw)[,c('condition', 'ourfilter')]), distfun='spearman')


filterCrit <- with(colData(scaRaw), pastFastqc=="PASS"& exonRate >0.3 & PercentToHuman>0.6 & nGeneOn> 4000)
sca <- subset(scaRaw,filterCrit)
sca <- sca[sample(which(freq(sca)>0), 6000),]

cdr2 <-colSums(assay(sca)>0)
qplot(x=cdr2, y=colData(sca)$nGeneOn) + xlab('New CDR') + ylab('Old CDR')

colData(sca)$cngeneson <- scale(cdr2)

scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=value))+geom_density() + facet_wrap(~symbolid, scale='free_y')


maits$expressionmat[1:10,1:10]
t(as.matrix(b.cells))[1:10,1:3]
dim(b.cells)

b.cells.sample
scaSample

flat[1:10]
b.cells.flat[1:10]

??Mast


