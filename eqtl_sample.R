# Load the clustered data.
load("../pbmc6k/data/seurat_33k_ensembl.Rda")

library("Matrix")
writeMM(pbmc33k.merged@raw.data, "~/Desktop/pbmc33k-raw-data.txt")
write(colnames(pbmc33k.merged@raw.data), file = "~/Desktop/pbmc33k-raw-data-cols.txt")
write(rownames(pbmc33k.merged@raw.data), file = "~/Desktop/pbmc33k-raw-data-rows.txt")

library("Matrix")
norm.33k <- readMM("~/Dropbox/Dropbox/pbmc33k/pbmc33k-norm-data.txt")
rownames(norm.33k) <- read.delim("~/Dropbox/Dropbox/pbmc33k/pbmc33k-norm-data-rows.txt", sep = "\n", header = F)[,1]
colnames(norm.33k) <- read.delim("~/Dropbox/Dropbox/pbmc33k/pbmc33k-norm-data-cols.txt", sep = "\n", header = F)[,1]
# Eventueel omzetten naar normale matrix
norm.33k <- as.matrix(norm.33k)

writeMM(pbmc33k.merged@data, "~/Desktop/pbmc33k-norm-data.txt")
write(colnames(pbmc33k.merged@data), file = "~/Desktop/pbmc33k-norm-data-cols.txt")
write(rownames(pbmc33k.merged@data), file = "~/Desktop/pbmc33k-norm-data-rows.txt")

sample.names <- c("1_LLDeep_0205","1_LLDeep_0209","1_LLDeep_0223","1_LLDeep_0229","1_LLDeep_0235","1_LLDeep_0240","1_LLDeep_0241","1_LLDeep_0253","1_LLDeep_0271","1_LLDeep_0288","1_LLDeep_0291","1_LLDeep_0292","1_LLDeep_0293","1_LLDeep_0298","1_LLDeep_0302","1_LLDeep_0341","1_LLDeep_0346","1_LLDeep_0358","1_LLDeep_0362","1_LLDeep_0364","1_LLDeep_0432","1_LLDeep_0624","1_LLDeep_0737","1_LLDeep_0755","1_LLDeep_0775","1_LLDeep_0783","1_LLDeep_0787","1_LLDeep_0862","1_LLDeep_0876","1_LLDeep_0880","1_LLDeep_0912","1_LLDeep_0938","1_LLDeep_1123","1_LLDeep_1156","1_LLDeep_1237")

nzmean <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else mean(x[!zvals])
}

nzmedian <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else median(x[!zvals])
}

nzmean.weighted <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else mean(x[!zvals]) * sum(!zvals) / length(x)
}

nzmedian.weighted <- function(x) {
  zvals <- x==0
  if (all(zvals)) 0 else median(x[!zvals]) * sum(!zvals) / length(x)
}


create.samples <- function(exp.data, sample.names) {
  exp.data <- exp.data[apply(exp.data, 1, function(x) !all(x==0)),]

  samples.matrix <- matrix(0L, nrow = nrow(exp.data), ncol = length(sample.names), dimnames = list(row.names(exp.data), sample.names))
  
  for (i in 1:length(sample.names)) {
    start <- (i * 500) - 499
    end <- i * 500
    samples.matrix[,i] <- rowMeans(exp.data[,start:end])
  }
  
  return(samples.matrix)
}

## Write the random samples to a file
samples.matrix <- create.samples(exp.data = as.matrix(pbmc33k.merged@data), sample.names = sample.names)
samples.matrix[1:10,1:5]
dim(samples.matrix)
write.table(samples.matrix, file = "traitfile-test.tsv", sep = "\t", quote = F)

## Boxplots and histogram of some samples 
hist(samples.matrix[,3], breaks = 100)
hist(samples.matrix[,3][samples.matrix[,3] != 0], breaks = 100)
boxplot(samples.matrix[,1:10], outline = F)

## Amount of zero values for 1 averaged sample
length(samples.matrix[,1])
sum(samples.matrix[,1] == 0)
sum(samples.matrix[,1] != 0)

## Histogram of 1 sample without genes with an average of 0
sample1 <- samples.matrix[,1]
sample1.not.null <- sample1[sample1 != 0]
hist(sample1.not.null, breaks = 100, ylim = c(0,200))
hist(samples.matrix[,3], breaks = 100, ylim = c(0,200))
## Compare to original
boxplot(sample1.not.null, samples.matrix[,1], outline = F)

sample.mean <- apply(pbmc33k.merged@data[,1:500], 1, mean)
sample.median <- apply(pbmc33k.merged@data[,1:500], 1, median)
sample.median.non.zero <- apply(pbmc33k.merged@data[,1:500],1,nzmedian)
sample.median.non.zero.weighted <- apply(pbmc33k.merged@data[,1:500],1,nzmedian.weighted)

## Calculate average of non-zero values random sample (n=500)
sample.mean.non.zero <- apply(pbmc33k.merged@data[,1:500],1,nzmean)
boxplot(sample.mean.non.zero)
boxplot(sample.mean.non.zero[sample.mean.non.zero != 0])


sample.mean.non.zero.weighted <- apply(pbmc33k.merged@data[,1:500],1,nzmean.weighted)
boxplot(sample.mean.non.zero.weighted)
boxplot(sample.mean[sample.mean != 0],
        sample.median,
        sample.mean.non.zero.weighted[sample.mean.non.zero.weighted != 0],
        sample.median.non.zero.weighted[sample.mean.non.zero.weighted != 0],
        sample.mean.non.zero[sample.mean.non.zero != 0],
        sample.median.non.zero[sample.median.non.zero != 0],
        outline = F,
        names = c("Mean", "Median", "Mean nz*frac", "Median nz*frac", "Mean nz", "Median nz"))

boxplot(sample.mean[sample.mean != 0],
        sample.median,
        sample.mean.non.zero.weighted[sample.mean.non.zero.weighted != 0],
        sample.median.non.zero.weighted[sample.mean.non.zero.weighted != 0],
        outline = F,
        names = c("Mean", "Median", "Mean nz*frac", "Median nz*frac"))

sample.sums <- apply(pbmc33k.merged@data[,1:500],1,sum)
boxplot(sample.sums[sample.sums!=0], outline=F)

## Check if non-zero are equal
length(sample1.not.null)
length(sample.mean.non.zero.weighted[sample.mean.non.zero.weighted != 0])
length(sample.mean.non.zero[sample.mean.non.zero != 0])

length()

sce <- newSCESet(countData=pbmc33k.merged@data[,1:500])
sce <- computeSumFactors(sce, positive=T)
summary(sizeFactors(sce))
sce <- normalize(sce)
exp <- norm_exprs(sce)
exp[1:10,1:10]

library(scran)



