pc1_2 <- read.delim("~/Documents/Pilot 3/eQTL//pc1_2.txt")

plot(pc1_2$Comp1, pc1_2$Comp2, xlab = "Component 1", ylab = "Component 2")

pc1_2[pc1_2$Comp1 < 0.1486,]
#
## Please fill in the outlier cut-off below
#
outlierCutoff <- 0
  
pdf("OutlierDetection.pdf")
palette(c("dodgerblue2", "red2"))
plot(pc1_2$Comp1, pc1_2$Comp2, xlab = "Component 1", ylab = "Component 2", col = as.numeric(pc1_2$Comp1 < outlierCutoff)+1, main = paste("Removed", sum(pc1_2$Comp1 < outlierCutoff), "expression outlier samples"))
abline(v=outlierCutoff, col = "red2")
dev.off()

## NB!: please check on which side of the outlierCutoff the outlier samples are:
# If outlier samples are on the left of the cutoff line run:
write.table(pc1_2$Sample[pc1_2$Comp1 >= outlierCutoff], file = "~/Documents/Pilot 3/eQTL/includedExpressionSamples.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# If outlier samples are on the right of the cutoff line run:
# write.table(pc1_2$Sample[pc1_2$Comp1 <= outlierCutoff], file = "includedExpressionSamples.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)