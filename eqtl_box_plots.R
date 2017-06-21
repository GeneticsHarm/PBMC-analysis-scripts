library(beeswarm)
library(ggplot2)

genotypes <- read.table("~/Documents/Pilot 3/eQTL/all_pairs_per_cell_type/eqtls_all_celltypes.genotypes.txt", check.names = F)
pbmc.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_all_log_norm.tsv", header = T, sep="\t", check.names = F)
b.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_b_log_norm.tsv", header = T, sep="\t", check.names = F)
b.scaled.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_b_scaled.tsv", header = T, sep="\t", check.names = F)
th.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_th_log_norm.tsv", header = T, sep="\t", check.names = F)
tc.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_tc_log_norm.tsv", header = T, sep="\t", check.names = F)
nk.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_nk_log_norm.tsv", header = T, sep="\t", check.names = F)
dend.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_dend_log_norm.tsv", header = T, sep="\t", check.names = F)
mono.expression <- read.table("~/Documents/R/Stage/data/pilot3_traitfile_mono_log_norm.tsv", header = T, sep="\t", check.names = F)

plot.eqtls <- function(genotypes, genotype, expression, gene, main ="") {
  genotypes <- genotypes[,colnames(genotypes) %in% colnames(expression)]
  genotypes <- droplevels(unlist(unname(genotypes[genotype,])))
  expression <- unlist(expression[gene,])
  #boxplot(expression~genotypes, outline = FALSE, main = paste0(gene, " (", genotype, ")"))
  boxplot(expression~genotypes, outline = FALSE, main = main)
  beeswarm(expression~genotypes, col = 4, pch = 16, add = TRUE)
}

genotypes <- genotypes[,colnames(genotypes) %in% colnames(pbmc.expression)]
genotypes[1:10,]

plot.eqtls(genotypes,"rs4821670", pbmc.expression, "LGALS2", main = "PBMC")
plot.eqtls(genotypes,"rs4821670", mono.expression, "LGALS2", main = "Mono")
plot.eqtls(genotypes,"rs4821670", th.expression, "LGALS2", main = "CD4+ T")
plot.eqtls(genotypes,"rs4821670", tc.expression, "LGALS2", main = "CD8+ T")
plot.eqtls(genotypes,"rs4821670", nk.expression, "LGALS2", main = "NK")
plot.eqtls(genotypes,"rs4821670", dend.expression, "LGALS2", main = "Dendritic")
plot.eqtls(genotypes,"rs4821670", b.expression, "LGALS2", main = "B")

plot.eqtls(genotypes,"rs1131017", pbmc.expression, "RPS26")
plot.eqtls(genotypes,"rs1131017", mono.expression, "RPS26")
plot.eqtls(genotypes,"rs1131017", th.expression, "RPS26")
plot.eqtls(genotypes,"rs1131017", tc.expression, "RPS26")
plot.eqtls(genotypes,"rs1131017", nk.expression, "RPS26")
plot.eqtls(genotypes,"rs1131017", dend.expression, "RPS26")

# B-cells
par(mfrow=c(2,3))
plot.eqtls(genotypes,"rs9332431", b.expression, "CHTF8", main = "B-cells")
plot.eqtls(genotypes,"rs9332431", pbmc.expression, "CHTF8", main = "PBMC")
plot.eqtls(genotypes,"rs9332431", th.expression, "CHTF8", main = "CD4+ T")
plot.eqtls(genotypes,"rs9332431", tc.expression, "CHTF8", main = "CD8+ T")
plot.eqtls(genotypes,"rs9332431", mono.expression, "CHTF8", main = "Monocytes")
plot.eqtls(genotypes,"rs9332431", nk.expression, "CHTF8", main = "NK")
plot.eqtls(genotypes,"rs9332431", dend.expression, "CHTF8", main = "Dendritic")
par(mfrow=c(1,1))


plot.eqtls(genotypes,"rs613854", b.expression, "RP5-1086K13.1")
plot.eqtls(genotypes,"rs613854", b.scaled.expression, "RP5-1086K13.1")


rs9332431 <- droplevels(unlist(unname(genotypes["rs9332431",])))
CHTF8 <- unlist(b.expression["CHTF8",])

boxplot(CHTF8~rs9332431, outline = FALSE)
rs9332431 <- rs9332431[names(pbmc.expression) %in% colnames(b.expression)]
data <- cbind.data.frame(CHTF8,rs9332431)
data

plot.all <- function(snp, gene) {
  snp.name <- snp
  snp <- droplevels(unlist(unname(genotypes[snp,])))
  
  snp = factor(snp, levels(snp)[c(3,2,1)])
  
  data <- data.frame(snp= snp, expression= unlist(pbmc.expression[gene,]), cell.type= "PBMC", color="black")
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(th.expression[gene,]), cell.type= "CD 4+ T", color="#153057"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(tc.expression[gene,]), cell.type= "CD 8+ T", color="#009ddb"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(nk.expression[gene,]), cell.type= "NK", color="#e64b50"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(mono.expression[gene,]), cell.type= "Monocytes", color="#edba1b"))
  data <- rbind.data.frame(data, data.frame(snp= snp[colnames(pbmc.expression) %in% colnames(b.expression)], expression= unlist(b.expression[gene,]), cell.type= "B", color="#71bc4b"))
  data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(dend.expression[gene,]), cell.type= "DC", color="#965ec8"))
  
  colors = c(rep("lightgrey", 3), rep("#153057",3), rep("#009ddb",3), rep("#e64b50",3), rep("#edba1b",3), rep("#71bc4b",3), rep("#965ec8",3))
  colors.strip = c("black","#153057","#009ddb","#e64b50","#edba1b","#71bc4b","#965ec8")
  

  ggplot(data, aes(x=snp, y=expression, group=snp)) +
    geom_boxplot(notch=F, color = "black", outlier.shape=NA, fill= colors, lwd=0.6, alpha=1) + 
    theme_minimal(base_family = "Helvetica Neue Light") +
    theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica Neue"),
          title = element_text(size = 20),
          axis.title.y = element_text(size = 12, family = "Helvetica Neue")) +
    #theme_bw(base_family = 'Helvetica') +
    #theme(panel.background = element_rect(fill = "#f6f6f7")) +
    facet_wrap(~cell.type, ncol = length(levels(data$cell.type)) ) +
    geom_point(position = position_jitter(width = 0.2), size = 0.8, color="black", alpha = 0.5) +
    ggtitle(paste0(gene, " / ", snp.name)) +
    ylab("Expression") +
    xlab("")
}

plot.all("rs4821670", "LGALS2")
plot.all("rs9332431", "CHTF8")
plot.all("rs9332431", "SNTB2")
plot.all("rs17184588", "PNRC2")
plot.all("rs1469335", "LILRA1")
plot.all("rs1645788", "LILRA3")
plot.all("rs6912538", "C6orf211")
plot.all("rs17184588", "PNRC2")
C6orf211

# Dendtritic cells
plot.all("rs28636077", "EIF5A")
plot.all("rs11123895", "RFX8")

# NK cells
plot.all("rs9469779", "RPS10")
plot.all("rs11658072", "EIF5A")
plot.all("rs1860596", "FAM184B")

# Tc's
plot.all("rs2301521", "SMDT1")
plot.all("rs28636077", "EIF5A")
plot.all("rs9368810", "RPS10")

# PBMC
plot.all("rs17184588", "PNRC2")


## Monique verzoek
plot.all("rs1131017", "RPS26")
plot.all("rs9469779", "RPS10")
plot.all("rs113166067", "KANSL1-AS1")
plot.all("rs11658072", "EIF5A")
plot.all("rs9271644", "HLA-DQA2")

plot.all("rs8141621", "LGALS2")

genotypes[,1:5]

snp <- droplevels(unlist(unname(genotypes["rs4821670",])))
expression <- unlist(mono.expression["LGALS2",])
pbmc <- unlist(pbmc.expression["LGALS2",])
data <- data.frame(snp= snp, expression= pbmc, cell.type= "PBMC")
data <- rbind.data.frame(data, data.frame(snp= snp, expression= expression, cell.type= "Monocytes"))
data <- rbind.data.frame(data, data.frame(snp= snp, expression= unlist(b.expression["LGALS2",]), cell.type= "B"))

data.frame(snp= snp, expression= unlist(b.expression["LGALS2",]), cell.type= "B")
snp
data

ggplot(data, aes(x=snp, y=expression, group=snp)) +
  geom_boxplot(notch=F, outlier.alpha = 0.3) + 
  theme_bw(base_family = 'Helvetica') +
  theme(panel.background = element_rect(fill = "#f6f6f7")) +
  facet_wrap(~cell.type)
beeswarm(CHTF8~rs9332431, col = 4, pch = 16, add = TRUE)


nk.expression[1:10,]

rs1131017 <- factor(rs1131017, levels(rs1131017)[c(2,1,3)])
RPS26 <- unlist(all.cells.expression["RPS26",])
boxplot(RPS26~rs1131017)
boxplot(RPS26~rs1131017, outline = FALSE, main = "RPS26 (rs1131017)")
beeswarm(RPS26~rs1131017, col = 4, pch = 16, add = TRUE)

## Second negative allele

pbmc.eqtls.centered <- read.table("~/Documents/Pilot 3/eQTL/all_pairs_per_cell_type/pbmc_eQTLProbesFDR0.05-ProbeLevel.txt", header = T, sep ="\t") 


strsplit()
strsplit(as.character(), split = '/')

check.false <- function(data) {
  data[sapply(strsplit(as.character(data$SNPType),"/"), `[`, 2) == data$AlleleAssessed & data$OverallZScore < 0,] 
}

th.eqtls.centered[sapply(strsplit(as.character(th.eqtls.centered$SNPType),"/"), `[`, 2) == th.eqtls.centered$AlleleAssessed & th.eqtls.centered$OverallZScore < 0,] 

pbmc.false <- check.false(pbmc.eqtls.centered)
pbmc.false[,c("SNPName", "HGNCName", "AlleleAssessed")]

plot.all("rs9469779", "RPS10")
plot.all("1:28237556", "THEMIS2")
sapply(strsplit(th.eqtls.centered$SNPTyp,"/"), `[`, 1)

genotypes["rs9469779",]

th.expression["RPS10","1_LLDeep_1123"]

rownames(genotypes)

genotypes["rs58220449",]

rownames(genotypes)

plot.all("5:64925256", "CENPK")


