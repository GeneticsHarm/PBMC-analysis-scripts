library(Seurat)

##
## Data loading
##

get.data <- function(pilot, lane) {
  data <- read.table(paste0("~/Documents/Pilot ",pilot,"/deAnonymizer/pilot",pilot,"lane",lane,"_results_deanonymize_chrALL.txt"), header = T, sep = "\t", fill=T)
  rownames(data) <- paste0(sub("(.+)\\-{1}\\d{1}", "\\1", data$cell_id), "_lane", lane)
  data$lane <- lane
  sample.names <- get.sample.names(data)
  data$assigned_sample.s. <- set.sample.name(data$Singlet_samp, sample.names)
  data.xy <- read.table(paste0("~/Documents/Pilot ",pilot,"/xy/lane_",lane,".txt"), header = T)
  data.xy <- data.xy[order(data.xy$cell),] # order
  data.exp <- Read10X(paste0("~/Documents/Pilot ",pilot,"/data/lane_",lane,"/"))
  data.seurat <- new("seurat", raw.data = data.exp)
  data.seurat <- Setup(data.seurat, min.cells = 0, min.genes = 0, project = paste0("pilot",pilot,".lane",lane), do.scale = F, do.center = F, names.field = 1, names.delim = "\\-")
  data.meta <- FetchData(data.seurat, c("nUMI", "nGene"))
  return(cbind(data, data.xy, data.meta))
}

get.sample.names <- function(data) {
  sample.names <- (c(as.character(data[data$outcome=="Singlet" & data$Singlet_samp == 0,]$assigned_sample.s.)[1],
                     as.character(data[data$outcome=="Singlet" & data$Singlet_samp == 1,]$assigned_sample.s.)[1],
                     as.character(data[data$outcome=="Singlet" & data$Singlet_samp == 2,]$assigned_sample.s.)[1],
                     as.character(data[data$outcome=="Singlet" & data$Singlet_samp == 3,]$assigned_sample.s.)[1],
                     as.character(data[data$outcome=="Singlet" & data$Singlet_samp == 4,]$assigned_sample.s.)[1],
                     as.character(data[data$outcome=="Singlet" & data$Singlet_samp == 5,]$assigned_sample.s.)[1]))
  return(sample.names)
}

set.sample.name <- function(id, names) names[id+1]

##
## Plot helper functions
##

source("./multiplot.R")

colors.x = c("#153057", "#009ddb", "#edba1b", "#71bc4b", "#e64b50", "#965ec8")

give.n <- function(x) {
  return(data.frame(y = unname(quantile(x, probs = 0.75, type = 3))-0.25, label = paste0("N=",length(x))))
}

boxplot.samples <- function(data, doublets=F, notch=T) {
  colors <- colors.x[!is.na(get.sample.names(data))]
  #if (!doublets) data <- subset(data, outcome=="Singlet")
  if (!doublets) data <- subset(data, llkDoublet.llkSinglet < 100)
  
  ggplot(data, aes(x=Singlet_samp, y=y, group=Singlet_samp)) +
    geom_boxplot(notch=notch, color = colors, outlier.alpha = 0.3) + 
    stat_summary(geom = 'point', fun.y=mean, shape = 23, size = 5, color = colors) +  
    theme_bw(base_family = 'Helvetica') +
    theme(panel.background = element_rect(fill = "#f6f6f7")) +
    guides(colour = guide_legend(override.aes = list(size=8))) +
    coord_cartesian(ylim = c(0, 12)) +
    stat_summary(fun.data = give.n, geom = "text", color = "black") +
    xlab("Samples") +
    ylab("Y Reads") +
    theme(plot.margin = unit(c(0.5,4.5,0.5,0.5), "cm"))
}

plot.scores <- function(data) {
  ggplot(data, aes(x=llkDoublet.llkSinglet, y=nGene, color = factor(Singlet_samp))) +
    geom_point(alpha = 1, size = 0.6) +
    theme_bw(base_size = 13, base_family = 'Helvetica') + 
    theme(legend.title=element_blank()) +
    scale_color_manual(labels=get.sample.names(data), values=colors.x) +
    guides(colour = guide_legend(override.aes = list(size=8))) +
    theme(panel.background = element_rect(fill = "#f6f6f7")) +
    geom_vline(xintercept = 0, color="red") +
    geom_vline(xintercept = 100, color="black") 
    #coord_cartesian(ylim = c(0, 3000), xlim = c(-250, 750))
}

plot.scores.dif <- function(data) {
  ggplot(data, aes(x=llkDoublet.llkSinglet, y=difSingletSamp, color = factor(Singlet_samp))) +
    geom_point(alpha = 1, size = 0.6) +
    theme_bw(base_size = 13, base_family = 'Helvetica') + 
    theme(legend.title=element_blank()) +
    scale_color_manual(labels=get.sample.names(data), values=colors.x) +
    guides(colour = guide_legend(override.aes = list(size=8))) +
    theme(panel.background = element_rect(fill = "#f6f6f7")) +
    geom_vline(xintercept = 0, color="red") +
    geom_vline(xintercept = 100, color="black") +
    geom_hline(yintercept = 2, color="black")
    #coord_cartesian(ylim = c(0, 3000), xlim = c(-250, 750))
}


plot.scores.ngenes <- function(data) {
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  ggplot(data, aes(x=llkDoublet.llkSinglet, y=llkDoublet, color = nGene)) +
    geom_point(alpha = 1, size = 0.6) +
    theme_bw(base_size = 13) + 
    scale_colour_gradientn(colours = myPalette(100)) +
    theme(panel.background = element_rect(fill = "#f6f6f7")) +
    geom_vline(xintercept = 0, color="red") +
    geom_vline(xintercept = 100, color="black") 
    #coord_cartesian(ylim = c(-1000, 0), xlim = c(-250, 750))
}

##
## Workspace
##

pilot4.lane1 <- get.data(pilot = 4, lane = 1)
pilot4.lane2 <- get.data(pilot = 4, lane = 2)
pilot4.lane3 <- get.data(pilot = 4, lane = 3)


pilot3.lane1 <- get.data(pilot = 3, lane = 1)
pilot3.lane2 <- get.data(pilot = 3, lane = 2)
pilot3.lane3 <- get.data(pilot = 3, lane = 3)
pilot3.lane4 <- get.data(pilot = 3, lane = 4)
pilot3.lane5 <- get.data(pilot = 3, lane = 5)
pilot3.lane6 <- get.data(pilot = 3, lane = 6)
pilot3.lane7 <- get.data(pilot = 3, lane = 7)
pilot3.lane8 <- get.data(pilot = 3, lane = 8)


pilot3.de.anonymize <- rbind(pilot3.lane1, pilot3.lane2, pilot3.lane3, pilot3.lane4, pilot3.lane5, pilot3.lane6, pilot3.lane7, pilot3.lane8)

dim(pilot3.de.anonymize)

divide.larger <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  return(max(a, b) / min(a, b))
}

pilot4.lane1$difSingletSamp <- mapply(divide.larger, pilot4.lane1$llkSingletSamp1, pilot4.lane1$llkSingletSamp2)
pilot4.lane1$difSingletSamp <- mapply(divide.larger, pilot4.lane1$llkSingletSamp1, pilot4.lane1$llkSingletSamp2)
pilot4.lane2$difSingletSamp <- mapply(divide.larger, pilot4.lane2$llkSingletSamp1, pilot4.lane2$llkSingletSamp2)
pilot4.lane1.gonl$difSingletSamp <- mapply(divide.larger, pilot4.lane1.gonl$llkSingletSamp1, pilot4.lane1.gonl$llkSingletSamp2)

pilot3.lane2$difSingletSamp <- mapply(divide.larger, pilot3.lane2$llkSingletSamp1, pilot3.lane2$llkSingletSamp2)
pilot3.lane6$difSingletSamp <- mapply(divide.larger, pilot3.lane6$llkSingletSamp1, pilot3.lane6$llkSingletSamp2)

plot.scores.dif(pilot4.lane1)
plot.scores.dif(pilot4.lane1.gonl)
plot.scores.dif(pilot4.lane2)
plot.scores.dif(pilot3.lane2)
plot.scores.dif(pilot3.lane6)

pilot4.lane1.gonl <- read.table("~/Documents/Pilot 4/deAnonymizer/pilot4lane1_gonl_results_deanonymize_chrALL.txt", header = T, sep = "\t", fill=T)
boxplot(pilot4.lane1.gonl$nSNPs_tested, pilot4.lane1[pilot4.lane1$outcome!="Singlet",]$nSNPs_tested)

multiplot(plot.scores(pilot4.lane1), boxplot.samples(pilot4.lane1), cols = 1)
multiplot(plot.scores(pilot4.lane2), boxplot.samples(pilot4.lane2), cols = 1)
multiplot(plot.scores(pilot4.lane3), boxplot.samples(pilot4.lane3, doublets = F), cols = 1)

multiplot(plot.scores(pilot3.lane1), boxplot.samples(pilot3.lane1, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane2), boxplot.samples(pilot3.lane2, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane3), boxplot.samples(pilot3.lane3, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane4), boxplot.samples(pilot3.lane4, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane5), boxplot.samples(pilot3.lane5, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane6), boxplot.samples(pilot3.lane6, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane7), boxplot.samples(pilot3.lane7, doublets = F), cols = 1)
multiplot(plot.scores(pilot3.lane8), boxplot.samples(pilot3.lane8, doublets = F), cols = 1)

boxplot(pilot3.lane2[pilot3.lane2$outcome=="Doublet",]$y)

pilot4.lane3.doubt <- pilot4.lane3[pilot4.lane3$llkDoublet.llkSinglet > 0 & pilot4.lane3$llkDoublet.llkSinglet < 100,]
pilot4.lane2.doubt <- pilot4.lane2[pilot4.lane2$llkDoublet.llkSinglet > 0 & pilot4.lane2$llkDoublet.llkSinglet < 100,]
pilot4.lane1.doubt <- pilot4.lane1[pilot4.lane1$llkDoublet.llkSinglet > 0 & pilot4.lane1$llkDoublet.llkSinglet < 100,]
pilot4.lane3.doubt$llkSingletSamp1/pilot4.lane3.doubt$llkSingletSamp2

barplot(table(c(pilot4.lane3.doubt$Singlet_samp, pilot4.lane2.doubt$Singlet_samp)))
barplot(table(c(pilot4.lane3.doubt$Doub_samp1, pilot4.lane3.doubt$Doub_samp2, pilot4.lane2.doubt$Doub_samp1, pilot4.lane2.doubt$Doub_samp2)))
barplot(table(pilot4.lane1.doubt$Singlet_samp))

pilot3.lane4.doubt <- pilot3.lane4[pilot3.lane4$llkDoublet.llkSinglet > 0 & pilot3.lane4$llkDoublet.llkSinglet < 100,]
pilot3.lane4.sure <- pilot3.lane4[pilot3.lane4$llkDoublet.llkSinglet > 100,]

boxplot(pilot3.lane4.doubt$nUMI / pilot3.lane4.doubt$nGene, pilot3.lane4.sure$nUMI/pilot3.lane4.sure$nGene,
        names = c("False positives", "True doublets"), ylab="UMIs / Gene", outline = F)

pilot3.lane4.doubt <- pilot3.lane4[pilot3.lane4$llkDoublet.llkSinglet > 0,]
nrow(pilot3.lane4.doubt) / nrow(pilot3.lane4)
nrow(pilot3.lane4.sure) / nrow(pilot3.lane4)

pilot3.lane1.doubt <- pilot3.lane1[pilot3.lane1$llkDoublet.llkSinglet > 0,]
pilot3.lane1.sure <- pilot3.lane1[pilot3.lane1$llkDoublet.llkSinglet > 100,]
nrow(pilot3.lane1.doubt) / nrow(pilot3.lane1)
nrow(pilot3.lane1.sure) / nrow(pilot3.lane1)

