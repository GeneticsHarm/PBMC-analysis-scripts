##
## Script for plotting networks graphs of the disease genes
## and an itterative test of the wilcoxon rank sum test
##
## Loaded data: ./data/pbmc33k.Rda
##

#load("./seurat_pbmc_6k.RData")
load("./data/pbmc33k.Rda")
load("./data/wilcox_tests.RData")

library(extrafont)
library(psych)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(scales)

create.network <- function(exp.data, genes) {
  # filter genes from expression data
  exp.data.genes <- as.matrix(exp.data[genes[genes %in% rownames(exp.data)],])
  # remove empty rows
  exp.data.genes <- exp.data.genes[apply(exp.data.genes[,-1], 1, function(x) !all(x==0)),]
  # calculate p-values
  exp.data.genes.cor.pvalue <- corr.test(t(as.matrix(exp.data.genes)), method = "pearson", adjust = "none")$p
  # create a network with an edge if p value < 0.05
  network.matrix <- ifelse(exp.data.genes.cor.pvalue < 0.05, 1, 0)
  
  node.sizes <- rescale(rowSums(exp.data.genes), to=c(0.2,20))
  
  # Draw the network
  network <- network(network.matrix, directed = F)
  plot <- ggnet2(network, label = T, color = "steelblue", size = node.sizes, label.size = 5, max_size = 20, alpha = 0.7, legend.position = "none")
  plot <- plot + theme(text=element_text(family="Wing"))
  
  print(plot)
  
  return(exp.data.genes.cor.pvalue)
}

wilcox.n.test <- function(x, y, rep=500) {
  p.values <- numeric(rep)
  for (i in 1:rep) p.values[i] <- wilcox.test(x, sample(y, length(x)))$p.value
  return(p.values)
}


TSNEPlot(pbmc33k)

#### Write data to matrix
total.data <- pbmc33k@data[apply(pbmc33k@data[,-1], 1, function(x) !all(x==0)),]
write.table(as.matrix(t(total.data)), file = "data/33k_matrix.txt", sep = "\t", quote = F)

#### Cell types
b.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("BCell_A", "BCell_B", "BCell_C"))]
b.cells <- b.cells[apply(b.cells[,-1], 1, function(x) !all(x==0)),]
b.cells[1:40,1:10]

b.cells.cor <- cor(t(as.matrix(b.cells)))
b.cells.r2 <- b.cells.cor ^ 2

b.cells.total.expression <- rowSums(as.matrix(b.cells))
b.cells.total.r2 <- rowSums(b.cells.r2)
b.cells.total.expression[b.cells.total.expression > 15000]

# Plot total expression against R2
plot(b.cells.total.r2~b.cells.total.expression)
# Plot R2 of a high expressed gene against total expression
plot(b.cells.r2["RPL10",]~b.cells.total.expression, ylim=c(0,0.3))
# Low expressed gene
plot(b.cells.r2[1,]~b.cells.total.expression, ylim=c(0,0.3))

b.cells.r2[1, b.cells.r2[1,] > 0.25]
b.cells["AL627309.1", b.cells["AL627309.1",] > 0]
b.cells["AC097713.4", b.cells["AC097713.4",] > 0]
b.cells["ANKRD31", b.cells["ANKRD31",] > 0]

not.it <- !(rownames(b.cells.cor) %in% diabetes.1.genes)


b.cells.not.dia <- b.cells[not.it,]
b.cells.not.dia.r2 <- b.cells.r2[not.it,not.it]

b.cells.dia <- b.cells[diabetes.1.genes[diabetes.1.genes %in% rownames(b.cells)],]
b.cells.dia.cor <- cor(t(as.matrix(b.cells.dia)))
b.cells.dia.r2 <- b.cells.dia.cor ^ 2

boxplot(b.cells.dia.cor, ylim=c(-0.05,0.06))

sd(b.cells.dia.r2)
sd(b.cells.not.dia.r2)

mean(as.vector(b.cells.dia))
mean(as.vector(b.cells.not.dia))

sd(b.cells.dia)
sd(b.cells.not.dia)

b.cells.dia.cor.upper <- b.cells.dia.cor[upper.tri(b.cells.dia.cor, diag = F)]
b.cells.not.dia.cor.upper <- b.cells.not.dia.cor[upper.tri(b.cells.not.dia.cor)]

### wilcoxon test against random samples (same size as disease gene pairs)
b.cells.not.dia.cor.upper.random <- sample(b.cells.not.dia.cor.upper, 1000000)
wilcox.test(b.cells.dia.cor.upper, b.cells.not.dia.cor.upper.random)
wilcox.test(b.cells.dia.cor.upper)

p.values <- wilcox.n.test(b.cells.dia.cor.upper, b.cells.not.dia.cor.upper)
mean(p.values)

pairwise.wilcox.test()

# wilcoxon test against all other samples
b.cells.wilcox <- wilcox.test(b.cells.dia.cor.upper, b.cells.not.dia.cor.upper)

t.h.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("CD4_Naive", "CD4_Memory_Antiviral", "CD4_Memory_TIGIT", "CD4_Memory"))]
t.h.cells <- t.h.cells[apply(t.h.cells[,-1], 1, function(x) !all(x==0)),]
t.h.cells.dia <-  as.matrix(t.h.cells[diabetes.1.genes[diabetes.1.genes %in% rownames(t.h.cells)],])
t.h.cells.dia.cor <- cor(t(as.matrix(t.h.cells.dia)))
t.h.cells.dia.cor.upper <- t.h.cells.dia.cor[upper.tri(t.h.cells.dia.cor)]

length(t.h.cells.dia.cor.upper)
length(b.cells.dia.cor.upper)
wilcox.test(t.h.cells.dia.cor.upper, b.cells.dia.cor.upper)


dim(t.h.cells)
t.c.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("CD8Effector_GZMH", "CD8Effector_GZMK", "CD8_Memory", "CD8_Naive"))]
dim(t.c.cells)
nk.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("NK", "NK_Prss57", "NK_XCL1"))]
dim(nk.cells)
monocytes <- pbmc33k@data[,WhichCells(pbmc33k, c("Mono_CD16+", "Mono_CD16+_C1qa", "Mono_Apobec3b", "Mono_Apobec3a", "Mono_CD14+", "Mono_CD14+_Inflam", "Mono_CD14+_Antiviral"))]
dim(monocytes)
megakaryos <- pbmc33k@data[,WhichCells(pbmc33k, "Megakaryocyte")]
dim(megakaryos)
dendritic.cells <- pbmc33k@data[,WhichCells(pbmc33k, c("DC_CD141+", "DC_CD1C+"))]
dim(dendritic.cells)

cor.test()

#### Disease genes

diabetes.1.genes <- c("AP4B1","BCAR1","C11orf21","C14orf64","C1QTNF6","C7orf71","CCDC101","CCR7","CD226","CDC37P1","CDK2","CLEC16A","CLECL1","CLN3","CSF3","CTLA4","CTSH","DOK6","EIF3C","EIF3CL","ENSG00000232530","ENSG00000241287","ENSG00000256582","ENSG00000257411","ENSG00000258511","ENSG00000261089","ENSG00000264198","ENSG00000266840","ERBB3","FKRP","GDF11","GSDMB","HLA.DMA","HLA.DQA1","HLA.DQB1","HLA.DQB1.AS1","HLA.DQB2","HLA.DRB5","HOTAIRM1","HOXA1","HOXA2","IKZF3","IL27","IL2RA","LIF","MIR320E","NPIPB6","NPIPB7","NPIPB8","NPIPB9","ORMDL3","PRKD2","PTPN2","RAB5B","RABEP2","RPS26","SBK1","SH2B1","SIRPG","SKAP2","STRN4","SULT1A1","SULT1A2","SUOX","TSPAN32","TYK2","ZPBP2")
celiac.genes <- c("AGER","ARHGAP31","BCL9L","C1orf106","CACNA1S","CCDC116","CCR1","CCR2","CCR3","CCRL2","CIITA","CPLX3","CSK","CXCR6","ELMO1","ENSG00000223552","ENSG00000236525","ENSG00000237058","ENSG00000238164","ENSG00000240310","ENSG00000241287","ENSG00000255422","ENSG00000259382","ENSG00000261606","ENSG00000263080","FLT1P1","GPR25","GPSM3","HLA.DMA","HLA.DQB1","HLA.DQB1.AS1","HLA.DRB5","ICOS","ICOSLG","IL18R1","IL1RL1","IRF4","KIAA1109","LMAN1L","LTF","MIR4513","MIR4772","MMEL1","PBX2","PLEK","PTPN2","RMI2","RNF5","THEMIS","TNFRSF14","TTC34","UBASH3A","UBE2L3","ULK3")
ibd.genes <- c("AAGAB","AAMP","ACO2","ACTR1A","ADO","AGAP5","ANKMY1","AP4B1","APEH","ARPC2","ATG16L1","BACH2","BSN","C10orf55","C17orf67","C20orf112","C5orf47","C7orf71","CACNA1S","CALM3","CAMK2G","CARD9","CAST","CCDC116","CCR5","CD226","CD40","CDH3","CEBPB","CELSR3","CISD1","CLN3","COMMD7","CPEB4","CRTC3","CSDC2","CTSW","CUL2","CXCL5","DACT3","DAG1","DENND1B","DESI1","DNLZ","DOCK7","DOK6","EDN3","EFCAB8","EFEMP2","EIF2S2P3","EIF3C","ENSG00000167807","ENSG00000170091","ENSG00000197146","ENSG00000202533","ENSG00000203999","ENSG00000213702","ENSG00000224276","ENSG00000224397","ENSG00000225331","ENSG00000225931","ENSG00000227598","ENSG00000228037","ENSG00000229299","ENSG00000229951","ENSG00000231858","ENSG00000231993","ENSG00000233077","ENSG00000234132","ENSG00000235582","ENSG00000236525","ENSG00000237604","ENSG00000238164","ENSG00000238280","ENSG00000242288","ENSG00000242687","ENSG00000243628","ENSG00000243696","ENSG00000247121","ENSG00000248636","ENSG00000249141","ENSG00000251572","ENSG00000254275","ENSG00000256428","ENSG00000256914","ENSG00000257582","ENSG00000261089","ENSG00000261338","ENSG00000261662","ENSG00000262820","ENSG00000266786","ENSG00000266840","ENSG00000267607","ENSG00000268810","ENSG00000269151","ENSG00000269292","ENSG00000269636","ERAP1","ERAP2","ERGIC1","FADS1","FADS2","FEN1","FIBP","FOXR1","FUT2","GNA12","GNG8","GPBAR1","GPR25","GSDMB","HLA.DRB5","HLA.DRB6","HOTAIRM1","HOXA1","HOXA2","ICAM1","ICAM4","ICAM5","ICOSLG","IL18R1","IL1RL1","IL23R","IL27","IL6ST","INPP5D","IP6K2","IPMK","IRF5","IRGM","ITGAL","ITIH4","ITM2BP1","KIF11","KIR2DL1","KIR3DL1","KIR3DL2","KSR1","LIME1","LNPEP","LRRK2","LY75","LY75.CD302","MAP3K8","MAPRE1","MEI1","MIR1908","MIR3939","MIR4772","MIR611","MMEL1","MRPL20","MRPL23","MST1","MUS81","MXRA8","MYRF","NDFIP1","NFATC1","NFKB1","NKX2.3","NOTCH2","NPIPB6","NPIPB7","NPIPB8","NPIPB9","NTN5","ORMDL3","P4HA2","PARK7","PEBP1P3","PF4","PFKFB4","PLAU","PLCH2","PLCL1","PMM1","POLR3H","PPP5C","PPP5D1","PRKAB1","PRKCB","PSMB8","PSMB9","PTGIR","PTGS2","PTPN2","PTPRC","RASGRP1","RBM22","REG4","RFTN2","RIPK2","RN7SL657P","RNASET2","RNF145","RNFT1","RPL9P29","RPS6KA2","SAG","SCNN1D","SEC1P","SKAP2","SLC22A5","SLC25A20","SLC26A6","SLC2A4RG","SLC35C2","SLC35D1","SLMO2","SMIM3","SMURF1","SNAPC4","SP140","SUFU","SULT1A2","SYNGR1","TEF","TM9SF4","TMBIM1","TMEM180","TMEM258","TNFRSF14","TNFSF8","TOB2","TSPAN14","TYK2","UBE2L3","UBQLN4","USP36","ZBTB38","ZBTB46","ZFP90","ZGLP1","ZGPAT","ZNF300P1","ZNF831","ZPBP2")
rheuma.genes <- c("ACSL6","AFF3","AHNAK2","ALS2CR12","ANKRD55","AP4B1","ARAP1","ARAP1-AS2","ARID5B","ATG5","BCL9L","BLK","CCDC157","CCR6","CD226 ","CD40","CD5","CFLAR","CTLA4","DDX6","ENSG00000228403","ENSG00000230955","ENSG00000237058","ENSG00000238164","ENSG00000254275","ENSG00000254774","ENSG00000255354","ENSG00000255518","ENSG00000256914","ENSG00000262211","ENSG00000264198","ENSG00000266840","ENSG00000267520","ENSG00000269468","ENSG00000269954","ERBB3","FADS1","FADS2","FAM167A","GPSM3","GSDMB","HLA-DQA1","HLA-DQA2","HLA-DRB1","HLA-DRB5","HLA-DRB6","HLA-DRB9","IFNAR1","IL6R","IL6ST","INPP5B","IRF5","JAZF1","KCTD20","MED1","MMEL1","MYRF","NOTCH4","ORMDL3","PADI4","PAM","PBX2","PDE2A","PHF19","PSMD5","RPS26","RTKN2","SHE","SLC35C2","SPRED2","SUOX","SYNGR1","TBC1D10A","TNFRSF14","TRAF1","TTC34","UBASH3A","UTS2","VPS37C","ZPBP2")

#### Diabetes
b.cells.diabetes.pvalues <- create.network(exp.data = b.cells, genes = diabetes.1.genes)
dim(b.cells.diabetes.pvalues)
t.h.cells.diabetes.pvalues <- create.network(exp.data = t.h.cells, genes = diabetes.1.genes)
dim(t.h.cells.diabetes.pvalues)
t.c.cells.diabetes.pvalues <- create.network(exp.data = t.c.cells, genes = diabetes.1.genes)
nk.cells.diabetes.pvalues <- create.network(exp.data = nk.cells, genes = diabetes.1.genes)
monocytes.diabetes.pvalues <- create.network(exp.data = monocytes, genes = diabetes.1.genes)
megakaryos.diabetes.pvalues <- create.network(exp.data = megakaryos, genes = diabetes.1.genes)
dendritic.cells.diabetes.pvalues <- create.network(exp.data = dendritic.cells, genes = diabetes.1.genes)
dim(dendritic.cells.diabetes.pvalues)

#### Rheuma
b.cells.rheuma.pvalues <- create.network(exp.data = b.cells, genes = rheuma.genes)
t.c.cells.rheuma.pvalues <- create.network(exp.data = t.c.cells, genes = rheuma.genes)
t.h.cells.rheuma.pvalues <- create.network(exp.data = t.h.cells, genes = rheuma.genes)

#### Inflammatory bowel disease
b.cells.ibd.pvalues <- create.network(exp.data = b.cells, genes = ibd.genes)
t.h.cells.ibd.pvalues <- create.network(exp.data = t.h.cells, genes = ibd.genes)
nk.cells.ibd.pavlues <- create.network(exp.data = nk.cells, genes = ibd.genes)

#### Celiac disease
b.cells.celiac.pvalues <- create.network(exp.data = b.cells, genes = celiac.genes)
t.h.cells.celiac.pvalues <- create.network(exp.data = t.h.cells, genes = celiac.genes)

