##
## Applies a Wilcoxon rank sum of the gene-pair correlations between diseases genes and non-diseases genes for every cell-type.
##
## Loaded data:
##     data/33k_correlation_spearman_celltype (created @ corr_spearman.R)
##     data/immune_disease_genes (created @ immune_diseases_genes.R)
##

load("data/33k_correlation_spearman_celltype")
load("data/immune_disease_genes")

## Helper function to execute the test between disease and non-disease gene pairs.
wilcox.diseases.genes <- function(pairs, disease.genes) {
  # Get the present disease and non disease genes pairs
  disease.gene.pairs <- pairs[pairs$gene1 %in% disease.genes & pairs$gene2 %in% disease.genes,]
  non.disease.gene.pairs <- pairs[!(pairs$gene1 %in% disease.genes & pairs$gene2 %in% disease.genes),]
  # Do the test
  test <- wilcox.test(disease.gene.pairs$rho^2, non.disease.gene.pairs$rho^2)
  # Add the present disease genes to wilcoxon.test object
  test$disease.genes <- unique(c(disease.gene.pairs$gene1, disease.gene.pairs$gene2))
  test$disease.gene.pairs <- disease.gene.pairs
  test$non.disease.gene.pairs <- non.disease.gene.pairs
  
  return(test)
}

compare.distribution <- function(pairs, disease.genes) {

  disease.gene.pairs <- pairs[pairs$gene1 %in% disease.genes & pairs$gene2 %in% disease.genes,]
  non.disease.gene.pairs <- pairs[!(pairs$gene1 %in% disease.genes & pairs$gene2 %in% disease.genes),]
  
  plot(density(non.disease.gene.pairs$rho^2))
  #lines(density(disease.gene.pairs$rho^2))
  
}

## Do the actual test
# Diabetes 1 genes
nk.dia.wilcox <- wilcox.diseases.genes(pairs = nk.cells.cor.pairs, disease.genes = diabetes.1.genes)
b.dia.wilcox <- wilcox.diseases.genes(pairs = b.cells.cor.pairs, disease.genes = diabetes.1.genes)
th.dia.wilcox <- wilcox.diseases.genes(pairs = t.h.cells.cor.pairs, disease.genes = diabetes.1.genes)
tc.dia.wilcox <- wilcox.diseases.genes(pairs = t.c.cells.cor.pairs, disease.genes = diabetes.1.genes)
mono.dia.wilcox <- wilcox.diseases.genes(pairs = monocytes.cor.pairs, disease.genes = diabetes.1.genes)
# Celiac disease
nk.celiac.wilcox <- wilcox.diseases.genes(pairs = nk.cells.cor.pairs, disease.genes = celiac.genes)
b.celiac.wilcox <- wilcox.diseases.genes(pairs = b.cells.cor.pairs, disease.genes = celiac.genes)
th.celiac.wilcox <- wilcox.diseases.genes(pairs = t.h.cells.cor.pairs, disease.genes = celiac.genes)
tc.celiac.wilcox <- wilcox.diseases.genes(pairs = t.c.cells.cor.pairs, disease.genes = celiac.genes)
mono.celiac.wilcox <- wilcox.diseases.genes(pairs = monocytes.cor.pairs, disease.genes = celiac.genes)
# Rheuma
nk.rheuma.wilcox <- wilcox.diseases.genes(pairs = nk.cells.cor.pairs, disease.genes = rheuma.genes)
b.rheuma.wilcox <- wilcox.diseases.genes(pairs = b.cells.cor.pairs, disease.genes = rheuma.genes)
th.rheuma.wilcox <- wilcox.diseases.genes(pairs = t.h.cells.cor.pairs, disease.genes = rheuma.genes)
tc.rheuma.wilcox <- wilcox.diseases.genes(pairs = t.c.cells.cor.pairs, disease.genes = rheuma.genes)
mono.rheuma.wilcox <- wilcox.diseases.genes(pairs = monocytes.cor.pairs, disease.genes = rheuma.genes)
# IBD
nk.ibd.wilcox <- wilcox.diseases.genes(pairs = nk.cells.cor.pairs, disease.genes = ibd.genes)
b.ibd.wilcox <- wilcox.diseases.genes(pairs = b.cells.cor.pairs, disease.genes = ibd.genes)
th.ibd.wilcox <- wilcox.diseases.genes(pairs = t.h.cells.cor.pairs, disease.genes = ibd.genes)
tc.ibd.wilcox <- wilcox.diseases.genes(pairs = t.c.cells.cor.pairs, disease.genes = ibd.genes)
mono.ibd.wilcox <- wilcox.diseases.genes(pairs = monocytes.cor.pairs, disease.genes = ibd.genes)

## Save the results
save(list = ls(pattern=".+wilcox"), file = "data/spearman_wilcox_test_frac_0.05.RData")
load("data/spearman_wilcox_test_frac_0.05.RData")

## Explore the results
nk.dia.wilcox$p.value
b.dia.wilcox$p.value
b.dia.wilcox$disease.gene.pairs
th.dia.wilcox$p.value
tc.dia.wilcox$p.value
mono.dia.wilcox$p.value
mono.dia.wilcox$disease.genes
mono.dia.wilcox$disease.gene.pairs
nk.dia.wilcox$disease.gene.pairs
b.dia.wilcox$disease.gene.pairs

nk.celiac.wilcox$p.value
b.celiac.wilcox$p.value
th.celiac.wilcox$p.value
tc.celiac.wilcox$p.value
tc.celiac.wilcox$disease.gene.pairs
mono.celiac.wilcox$p.value
mono.celiac.wilcox$disease.gene.pairs
mono.celiac.wilcox$disease.genes

nk.rheuma.wilcox$p.value
b.rheuma.wilcox$p.value
th.rheuma.wilcox$p.value
tc.rheuma.wilcox$p.value
mono.rheuma.wilcox$p.value
nk.rheuma.wilcox$disease.genes
nk.rheuma.wilcox$disease.gene.pairs

nk.ibd.wilcox$p.value
b.ibd.wilcox$p.value
th.ibd.wilcox$p.value
tc.ibd.wilcox$p.value
mono.ibd.wilcox$p.value

nk.ibd.wilcox$disease.genes
nk.ibd.wilcox$disease.gene.pairs
mono.ibd.wilcox$disease.gene.pairs[1:10,]
mono.dia.wilcox$disease.gene.pairs[1:10,]

#### WORKSPACE

b.cells.cor.pairs

sqrt(7263766*2)
3812*3811/2

compare.distribution(b.cells.cor.pairs, diabetes.1.genes)
h = hist(nk.cells.cor.pairs$rho)
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE)
# Add a Normal Curve (Thanks to Peter Dalgaard)
x <- nk.cells.cor.pairs$rho

h<-hist(x, breaks=1000) 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
