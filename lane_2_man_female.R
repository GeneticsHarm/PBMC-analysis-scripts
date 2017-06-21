load("data/pilot_3_clustered.Rda")
FeaturePlot(combined.seurat.subset, "XIST")

sum(combined.seurat.subset@raw.data[c("XIST", "FTX", "JPX"),] == 0)/length(combined.seurat.subset@raw.data[c("XIST", "FTX", "JPX"),])
combined.seurat.subset@raw.data[1:10,1:10]

hist(combined.seurat.subset@raw.data["XIST",], bre)

load("../")

lane1.xy <- read.table("~/Documents/Pilot 3/xy/lane_1.txt", header = T)
lane1.data <- Read10X("~/Documents/Pilot 3/data/lane_1/")

lane3.data <- Read10X("~/Documents/Pilot 3/data/lane_3/")
lane3.xy <- read.table("~/Documents/Pilot 3/xy/lane_3.txt", header = T)
lane3.ids <- read.table("~/Documents/Pilot 3/deAnonymizer/p3l3.txt", header = T)

lane3.xy[,"cell"] <- sub("([ATCG]+)-\\d{1}", "\\1", lane3.xy[,"cell"]) # trim the barcode
lane3.xy <- lane3.xy[order(lane3.xy$cell),] # order
rownames(lane3.ids) <- sub("([ATCG]+)-\\d{1}", "\\1", rownames(lane3.ids)) # trim the barcode

boxplot(lane3.xy$y~lane3.ids$samp, ylab = "Y reads", xlab = "Sample")

lane3.doublets <- lane3.ids$pos_doublet == 1
boxplot(lane3.xy$y[!lane3.doublets]~lane3.ids$samp[!lane3.doublets], ylab = "Y reads", xlab = "Sample")

length(lane3.ids[lane3.ids$samp == 2,]$samp)


lane2.xy <- read.table("~/Documents/Pilot 3/xy/lane_2.txt", header = T)
lane2.data <- Read10X("~/Documents/Pilot 3/data/lane_2/")
lane2.ids <- read.table("~/Documents/Pilot 3/deAnonymizer/p3l2.txt", header = T)
lane2.ids[1:10,]
lane2.xy[,"cell"] <- sub("([ATCG]+)-\\d{1}", "\\1", lane2.xy[,"cell"]) # trim the barcode
lane2.xy <- lane2.xy[order(lane2.xy$cell),] # order
rownames(lane2.ids) <- sub("([ATCG]+)-\\d{1}", "\\1", rownames(lane2.ids)) # trim the barcode
sum(rownames(lane2.ids) == lane2.xy[,"cell"])

boxplot(lane2.xy$y[!lane2.doublets]~lane2.ids$sample[!lane2.doublets], ylab = "Y reads", xlab = "Sample")
lane2.doublets <- lane2.ids$pos_doublet == 1
sum(lane2.doublets)

lane1.xy[,"cell"] <- sub("([ATCG]+)-\\d{1}", "\\1", lane1.xy[,"cell"]) # trim the barcode
lane1.xy <- lane1.xy[order(lane1.xy$cell),] # order
sum(lane1.xy$cell == colnames(lane1.data)) # check if barcodes are in equal order 

lane1.data[1:10,1:3]
lane1.data["XIST",]
plot(lane1.xy$y, lane1.data["XIST",], pch=16, col = "black")
dotchart(lane1.xy$y, lane1.data["XIST",])

ggplot(lane1.xy, aes(x = y , y = lane1.data["XIST",])) + geom_point(aes(size = y))

female.genes <- lane1.data["FTX",] + lane1.data["XIST",] + lane1.data["JPX",]

### base graphics ###
plot(lane1.xy$y, female.genes, pch = 16, cex = .9)
lane1.xy[1:10,]

boxplot()
plot(lane1.xy$y, lane1.xy$x/lane1.xy$total, pch = 16)
boxplot(lane1.xy$x/lane1.xy$total~lane1.xy$y, xlab="Aantal y reads", ylab = "X reads / total reads")
boxplot(lane1.xy$x/(lane1.xy$x+lane1.xy$y)~lane1.xy$y, xlab="Aantal y reads", ylab = "X reads / total sex chrom reads")
boxplot(lane1.xy$total~lane1.xy$y, xlab="Aantal y reads", ylab = "total reads")
plot(lane1.xy$y, lane1.xy$total)

