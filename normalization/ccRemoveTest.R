##
## Script used to regress out cell-cycle effects on scRNA-seq data
##


#install.packages("~/Downloads/ccRemover/ccRemover_0.3.tar.gz", repos=NULL)
library(ccRemover)

memory.limit()

# Load data.
pbmc.data <- read.table("~/Documents/Stage/Matrices/count-matrix-6k-v2.tsv", row.names=1, header=T)
# Remove all genes with zero counts.
pbmc.data <- pbmc.data[apply(pbmc.data[,-1], 1, function(x) !all(x==0)),]
# Normalize
pbmc.data.norm <- log2(pbmc.data + 1) 
# Scale to mean of 0 for each gene
pbmc.data.norm.scaled <- t(scale(t(pbmc.data.norm), scale=FALSE, center=TRUE))
# Check which genes are cell cycle genes
load("cell_cycle_gene_indexer/HScc_genes.rda")
if.cc <- rownames(pbmc.data.norm.scaled) %in% human_cell_cycle_genes[,1]

head(rownames(pbmc.data.norm.scaled))

sum(if.cc, na.rm=TRUE)

# Create a list with data vector
pbmc.data.cc <- list(x=pbmc.data.norm.scaled, if.cc=if.cc)
# Run ccRemover
pmbc.data.corrected <- ccRemover.main(pbmc.data.cc, nboot=25)

??ccRemover
?"Memory-limits"
