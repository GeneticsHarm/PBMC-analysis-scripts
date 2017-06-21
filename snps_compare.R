th.eqtls.log.norm <- read.table("~/Documents/Pilot 3/eQTL/th-cells/eQTLProbesFDR0.05-ProbeLevel_log_norm.txt", header = T, sep = "\t") 
th.eqtls.centered <- read.table("~/Documents/Pilot 3/eQTL/th-cells/eQTLProbesFDR0.05-ProbeLevel_scaled_centered.txt", header = T, sep ="\t") 
eqtls <- read.table("~/Documents/Pilot 3/eQTL/eQTLProbesFDR0.05-ProbeLevel.txt", header = T, sep = "\t")

hist(th.eqtls.centered$SNPChrPos - th.eqtls.centered$ProbeCenterChrPos, breaks = 100)
hist(th.eqtls.log.norm$SNPChrPos - th.eqtls.log.norm$ProbeCenterChrPos, breaks = 100)

##
## Compare pilot 3 snps with eQTL gen snps
##

snps.ll <- readLines("~/Documents/Pilot 3/eQTL/SNPs_LL_egcut.txt")
snps.pilot3 <- readLines("~/Documents/Pilot 3/eQTL/SNPs_pilot3.txt")
snps.pilot3.harmonized <- readLines("~/Documents/Pilot 3/eQTL/SNPs_pilot3_harmonized.txt")

snps.intersection <- snps.ll %in% snps.pilot3.harmonized
sum(snps.intersection)

##
## DeepSAGE vs. T-h centered
##

eqtls.deepsage <- read.table("~/Documents/Pilot 3/eQTL/DeepSAGE_results.csv", sep = ";", header = T, dec = ",")
#eqtls.th.deepsage <- read.table("~/Documents/Pilot 3/eQTL/DeepSAGE_eQTLsFDR-ProbeLevel.txt", header = T, sep = "\t")
eqtls.all.centered.deepsage <- read.table("~/Documents/Pilot 3/eQTL/DeepSAGE_all_cells_centered_eQTLsFDR-ProbeLevel.txt", header = T, sep = "\t")
#eqtls.all.deepsage <- read.table("~/Documents/Pilot 3/eQTL/DeepSAGE_all_cells_eQTLsFDR-ProbeLevel.txt", header = T, sep = "\t")

eqtls.deepsage$combine = paste0(eqtls.deepsage$SNPName, eqtls.deepsage$HGNCName)
eqtls.all.centered.deepsage$combine = paste0(eqtls.all.centered.deepsage$SNPName, eqtls.all.centered.deepsage$HGNCName)
#eqtls.all.deepsage$combine = as.character(interaction(eqtls.all.deepsage$SNPName, eqtls.all.deepsage$ProbeName))

#eqtls.deepsage <- eqtls.deepsage[eqtls.deepsage$combine %in% eqtls.th.deepsage$combine,]
eqtls.matched.deepsage <- eqtls.deepsage[match(eqtls.all.centered.deepsage$SNPName,eqtls.deepsage$SNPName, nomatch = 0),]
eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.centered.deepsage$AlleleAssessed,]$OverallZScore <- eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.centered.deepsage$AlleleAssessed,]$OverallZScore * -1

eqtls.matched.deepsage$SNP.probe.dist <- eqtls.matched.deepsage$SNPChrPos - eqtls.matched.deepsage$ProbeCenterChrPos
eqtls.th.centered.deepsage$SNP.probe.dist <- eqtls.th.centered.deepsage$SNPChrPos - eqtls.th.centered.deepsage$ProbeCenterChrPos

###
### Z-scores 
###
all.cells.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_all_log_norm.tsv", header = T, sep="\t")

exp.sums <- rowMeans(all.cells.expression)
exp.sums <- exp.sums[match(eqtls.all.centered.deepsage$HGNCName, names(exp.sums), nomatch = 0)]

plot(eqtls.matched.deepsage$OverallZScore[1:15], eqtls.all.centered.deepsage$OverallZScore[1:15], type = "p", 
     pch=21, bg=rgb(0,0,0,0.4), cex=exp.sums * 10, xlab = "DeepSAGE effect size (Z-scores)", ylab = "scRNA-seq effect size (Z-scores)", lwd=.8, main = "scRNA-seq (PBMC) vs. DeepSAGE (whole blood)", cex.lab = 1.4, cex.main=1.5, xlim = c(-10,10))
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)

layout(matrix(1:2,ncol=2),widths=c(2,2),heights=c(2,2),respect=FALSE)

par(mar = c(1,1,1,0), oma=c(4, 3, 2, 2), las=1)

plot(eqtls.matched.deepsage$OverallZScore, eqtls.all.centered.deepsage$OverallZScore, type = "p", 
     pch=21,
     col=c(rep("red", 15),rep("black",length(eqtls.all.centered.deepsage$OverallZScore))),
     bg=c(rep(rgb(1,0,0,1),15), rep(rgb(0,0,0,0.4), length(eqtls.all.centered.deepsage$OverallZScore))),
     cex=exp.sums + 0.5,
     xlab = "",
     ylab = "",
     lwd=.8,cex.lab = 1.4,
     xlim = c(-10,10))

abline(a=0, b=0, h=0,lty=1)
abline(a=0, b=0, v=0,lty=1)
abline(a=0, b=0, h=3.43,lty=3)
abline(a=0, b=0, h=-3.43,lty=3)

plot(eqtls.matched.deepsage$OverallZScore[1:15], eqtls.all.centered.deepsage$OverallZScore[1:15], type = "p", 
     pch=21, bg=rgb(0,0,0,0.4), cex=exp.sums * 10,
     xlab = "",
     ylab = "", lwd=.8,
     cex.lab = 1.4, cex.main=1.5, xlim = c(-10,10),
     yaxt='n', ann=FALSE)
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)

mtext("DeepSAGE effect size (Z-scores)", 1, 2, outer=TRUE, cex = 1)
mtext("scRNA-seq effect size (Z-scores)", 2, 2, outer=TRUE, las=0, cex =1)
mtext("A", 3, 0, adj=0, outer=TRUE, cex = 1.5)
mtext("B", 3, 0, outer=TRUE, cex = 1.5)


library(plotrix)
gap.plot(eqtls.matched.deepsage$OverallZScore, eqtls.th.centered.deepsage$OverallZScore,
         gap = c(-3.5,3.5), gap.axis = "x", bgcol = "white", xlim = c(-10,10),
         pch=21, bg=c(rep(rgb(1,0,0,1),19), rep(rgb(0,0,0,0.4),nrow(eqtls.matched.deepsage))), cex=exp.sums / 2 + 0.5, xlab = "DeepSAGE Z-scores", ylab = "Th cells Z-scores", lwd=.8, main = "CD4+ T cells vs. DeepSAGE")
abline(a=0, b=0, h=0,lty=3)
#abline(a=0, b=0, v=0,lty=3)

##
## Hist
##
hist(eqtls.matched.deepsage$SNP.probe.dist, breaks = 100)
hist(eqtls.th.centered.deepsage$SNP.probe.dist, breaks = 100)

eqtls.deepsage <- eqtls.deepsage[eqtls.deepsage$combine %in% eqtls.all.deepsage$combine,]
eqtls.matched.deepsage <- eqtls.deepsage[match(eqtls.all.deepsage$combine,eqtls.deepsage$combine, nomatch = 0),]
#eqtls.matched.deepsage <- eqtls.matched.deepsage[1:286,]
eqtls.matched.deepsage$combine == eqtls.all.deepsage$combine
eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.deepsage$AlleleAssessed,]$OverallZScore <- eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.deepsage$AlleleAssessed,]$OverallZScore * -1
plot(eqtls.matched.deepsage$OverallZScore[1:38], eqtls.all.deepsage$OverallZScore[1:38], pch=19, cex=0.5, xlab = "DeepSAGE Z-Scores", ylab = "All cells Z-Scores")
eqtls.all.deepsage[,1:11]

# Change direction when other allele is assessed
eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.th.deepsage$AlleleAssessed,]$OverallZScore <- eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.th.deepsage$AlleleAssessed,]$OverallZScore * -1

plot(eqtls.matched.deepsage$OverallZScore[1:28], eqtls.th.deepsage$OverallZScore[1:28], pch=19, cex=0.5, xlab = "DeepSAGE Z-Scores", ylab = "Th cells Z-Scores")
plot(eqtls.matched.deepsage$OverallZScore, eqtls.th.deepsage$OverallZScore, pch=19, cex=0.5, xlab = "DeepSAGE Z-Scores", ylab = "Th cells Z-Scores")


eqtls.matched.deepsage[1:10,1:11]
eqtls.th.deepsage[1:10,1:11]



length(snps.pilot3)
snps.ll[snps.intersection][1:500]
snps.ll[!snps.intersection][1:500]
snps.ll[1:200]
snps.pilot3[1:200]

length(unique(eqtls$SNPName))
eqtls.snps.in.pilot3 <- eqtls$SNPName %in% snps.pilot3 

length(unique(eqtls$SNPName))
sum(eqtls.snps.in.pilot3)
eqtls$SNPName[!eqtls.snps.in.pilot3]

eqtls[1:10,]

dim(eqtls)
dim(th.eqtls.centered)
dim(th.eqtls.log.norm)

sum(th.eqtls.centered$SNPName %in% eqtls$SNPName)
sum(th.eqtls.log.norm$SNPName %in% eqtls$SNPName)
sum(th.eqtls.log.norm$SNPName %in% th.eqtls.centered$SNPName)
th.eqtls.centered[1,1:10]
th.eqtls.log.norm[th.eqtls.log.norm$SNPName %in% eqtls$SNPName,1:11]
eqtls[eqtls$SNPName %in% th.eqtls.log.norm$SNPName,1:11]

##
## RNA-seq eQTLS
##
eqtls.rnaseq <- read.table("~/Documents/Pilot 3/eQTL/gene_level_eQTLs_independent_effects_interactions.txt", header = T, sep = "\t")
#eqtls.sc.rnaseq <- read.table("~/Documents/Pilot 3/eQTL/RNA-seq_all_cells_eQTLsFDR-ProbeLevel.txt", header = T, sep = "\t")
eqtls.sc.rnaseq <- read.table("~/Documents/Pilot 3/eQTL/all_centered_RNA-seq_eQTLsFDR-ProbeLevel.txt", header = T, sep = "\t")

eqtls.rnaseq$combine = paste(eqtls.rnaseq$SNP, eqtls.rnaseq$Gene)
eqtls.sc.rnaseq$combine = paste(eqtls.sc.rnaseq$SNPName, eqtls.sc.rnaseq$ProbeName)

eqtls.rnaseq$combine[1:10]
eqtls.sc.rnaseq$combine[1:10]

eqtls.rnaseq.matched <- eqtls.rnaseq[match(eqtls.sc.rnaseq$combine,eqtls.rnaseq$combine, nomatch = 0),]

eqtls.rnaseq.matched$Assesed.Allele == eqtls.sc.rnaseq$AlleleAssessed

eqtls.rnaseq.matched[eqtls.rnaseq.matched$Assesed.Allele != eqtls.sc.rnaseq$AlleleAssessed,]$Z.score <- eqtls.rnaseq.matched[eqtls.rnaseq.matched$Assesed.Allele != eqtls.sc.rnaseq$AlleleAssessed,]$Z.score * -1
plot(eqtls.rnaseq.matched$Z.score[1:48], eqtls.sc.rnaseq$OverallZScore[1:48])

all.cells.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_all_log_norm.tsv", header = T, sep="\t")
all.cells.expression[1:10,1:10]
genotypes <- read.table("~/Documents/Pilot 3/eQTL/genotypes.genotypes.txt")

genotypes <- genotypes[,colnames(genotypes) %in% colnames(all.cells.expression)]


rs1131017 <- unlist(unname(genotypes[1,]))
rs1131017 <- factor(rs1131017, levels(rs1131017)[c(2,1,3)])
RPS26 <- unlist(all.cells.expression["RPS26",])

boxplot(RPS26~rs1131017)

boxplot(RPS26~rs1131017, outline = FALSE, main = "RPS26 (rs1131017)")
beeswarm(RPS26~rs1131017, col = 4, pch = 16, add = TRUE)

exp.sums <- rowMeans(all.cells.expression)
exp.sums <- exp.sums[match(eqtls.rnaseq.matched$Gene.name, names(exp.sums), nomatch = 0)]
exp.sums[1:20]
is.30 <- (abs(eqtls.rnaseq.matched$Z.score) > 30)[1:90]
sum(is.30)


plot(eqtls.rnaseq.matched$Z.score[1:90][is.30], eqtls.sc.rnaseq$OverallZScore[1:90][is.30], type = "p", 
     pch=21, bg=rgb(0,0,0,0.4), cex=exp.sums[1:90][is.30] + 0.3, xlab = "RNA-seq Z-scores", ylab = "scRNA-seq Z-scores", lwd=.8, main = "All cell-types vs. BIOS RNA-seq")
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)

all.cells.expression[1:10,1:10]

eqtls.sc.rnaseq[1:90,][is.30,1:11]
eqtls.rnaseq.matched[1:90,][is.30,1:11]

## Th-cells centered vs RNA-seq top effects
eqtls.rnaseq <- read.table("~/Documents/Pilot 3/eQTL/gene_level_eQTLs_independent_effects_interactions.txt", header = T, sep = "\t")
eqtls.th.centered.rnaseq <- read.table("~/Documents/Pilot 3/eQTL/th_centered_RNA-seq_eQTLsFDR-ProbeLevel.txt", header = T, sep = "\t")
eqtls.th.centered.rnaseq$combine = paste(eqtls.th.centered.rnaseq$SNPName, eqtls.th.centered.rnaseq$ProbeName)

eqtls.level1.rnaseq <- eqtls.rnaseq[eqtls.rnaseq$eQTLLevelNumber == 1,]
eqtls.level1.rnaseq$combine = paste(eqtls.level1.rnaseq$SNP, eqtls.level1.rnaseq$Gene.name)
eqtls.th.centered.rnaseq <- eqtls.th.centered.rnaseq[eqtls.th.centered.rnaseq$combine %in% eqtls.level1.rnaseq$combine,]

eqtls.rnaseq.matched <- eqtls.level1.rnaseq[match(eqtls.th.centered.rnaseq$combine, eqtls.level1.rnaseq$combine, nomatch = 0),]



eqtls.rnaseq.matched[1:10,1:10]
eqtls.th.centered.rnaseq[1:10,1:10]

eqtls.rnaseq.matched[eqtls.rnaseq.matched$Assesed.Allele != eqtls.th.centered.rnaseq$AlleleAssessed,]$Z.score <- eqtls.rnaseq.matched[eqtls.rnaseq.matched$Assesed.Allele != eqtls.th.centered.rnaseq$AlleleAssessed,]$Z.score * -1

## Get gene expression
#th.cells.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_th_scaled_centered.tsv", header = T, sep="\t")
th.cells.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_th.tsv", header = T, sep="\t")

exp.sums <- rowMeans(th.cells.expression)
exp.sums <- exp.sums[match(eqtls.rnaseq.matched$Gene.name, names(exp.sums), nomatch = 0)]

exp.sums[1:10]

plot(eqtls.rnaseq.matched$Z.score[1:43], eqtls.th.centered.rnaseq$OverallZScore[1:43], type = "p", 
     pch=21, bg=rgb(0,0,0,0.4), cex=exp.sums[1:43] + 0.3, xlab = "RNA-seq Z-scores", ylab = "Th cells Z-scores", lwd=.8, main = "CD4+ T cells vs. BIOS RNA-seq")
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)

plot(eqtls.rnaseq.matched$Z.score, eqtls.th.centered.rnaseq$OverallZScore, type = "p",
     pch=21, bg=rgb(0,0,0,0.3), cex=exp.sums + 0.1, xlab = "RNA-seq Z-scores", ylab = "Th cells Z-scores", lwd=.5, main = "CD4+ T cells vs. BIOS RNA-seq")
abline(a=0, b=0, h=0,lty=3)
abline(a=0, b=0, v=0,lty=3)



hist(eqtls.rnaseq.matched$SNP.Chr.Position - eqtls.rnaseq.matched$Gene.Chr.position)
hist(eqtls.th.centered.rnaseq$SNPChrPos - eqtls.th.centered.rnaseq$ProbeCenterChrPos, breaks = 200, xlim = c(-250000,250000))
eqtls.th.centered.rnaseq$distance <- eqtls.th.centered.rnaseq$SNPChrPos - eqtls.th.centered.rnaseq$ProbeCenterChrPos
hist(eqtls.th.centered.rnaseq$distance[1:100], breaks = 100)

#eqtls.matched.deepsage <- eqtls.matched.deepsage[1:286,]
eqtls.matched.deepsage$combine == eqtls.all.deepsage$combine
eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.deepsage$AlleleAssessed,]$OverallZScore <- eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.all.deepsage$AlleleAssessed,]$OverallZScore * -1
plot(eqtls.matched.deepsage$OverallZScore[1:38], eqtls.all.deepsage$OverallZScore[1:38], pch=19, cex=0.5, xlab = "DeepSAGE Z-Scores", ylab = "All cells Z-Scores")
eqtls.all.deepsage[,1:11]

# Change direction when other allele is assessed
eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.th.deepsage$AlleleAssessed,]$OverallZScore <- eqtls.matched.deepsage[eqtls.matched.deepsage$AlleleAssessed != eqtls.th.deepsage$AlleleAssessed,]$OverallZScore * -1


##
##
eqtls.rnaseq[eqtls.rnaseq$Gene.Chr.position == 32627635,]






