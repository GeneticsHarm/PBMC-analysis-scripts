lane1.data <- Read10X("~/Documents/Pilot 3/data/lane_1/")
colnames(lane1.data) <- paste0(colnames(lane1.data), "_lane1")
lane2.data <- Read10X("~/Documents/Pilot 3/data/lane_2/")
colnames(lane2.data) <- paste0(colnames(lane2.data), "_lane2")
lane3.data <- Read10X("~/Documents/Pilot 3/data/lane_3/")
colnames(lane3.data) <- paste0(colnames(lane3.data), "_lane3")
lane4.data <- Read10X("~/Documents/Pilot 3/data/lane_4/")
colnames(lane4.data) <- paste0(colnames(lane4.data), "_lane4")
lane5.data <- Read10X("~/Documents/Pilot 3/data/lane_5/")
colnames(lane5.data) <- paste0(colnames(lane5.data), "_lane5")
lane6.data <- Read10X("~/Documents/Pilot 3/data/lane_6/")
colnames(lane6.data) <- paste0(colnames(lane6.data), "_lane6")
lane7.data <- Read10X("~/Documents/Pilot 3/data/lane_7/")
colnames(lane7.data) <- paste0(colnames(lane7.data), "_lane7")
lane8.data <- Read10X("~/Documents/Pilot 3/data/lane_8/")
colnames(lane8.data) <- paste0(colnames(lane8.data), "_lane8")

load(file = "./data/pilot_3.Rda")

pilot3.meta <- FetchData(combined.seurat, c("nUMI", "nGene", "percent.mito", "orig.ident"))

combined.mad <- mad(pilot3.meta$nGene)
nUMI.log.norm <- colSums(combined.seurat@data)
nUMI.log.norm.mad <- mad(nUMI.log.norm)
nUMI.cutoff <- nUMI.log.norm.mad * 3 + median(nUMI.log.norm)

nGene.mad <- mad(pilot3.meta$nGene)
nGene.cutoff <- nGene.mad * 3 + median(pilot3.meta$nGene)
nGene.sd <- sd(pilot3.meta$nGene)
nGene.high.cutoff.sd <-  median(pilot3.meta$nGene) + nGene.sd * 3
nGene.low.cutoff.sd <- median(pilot3.meta$nGene) - nGene.sd * 3

## Gene/reads per cell plot
ggplot(pilot3.meta, aes(nUMI.log.norm, nGene)) +
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Number of genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 2500, color="red")  + 
  geom_hline(yintercept = 500, color="red") +
  geom_vline(xintercept = nUMI.cutoff, color="green") + 
  geom_hline(yintercept = nGene.cutoff, color="green") +
  geom_hline(yintercept = nGene.high.cutoff.sd, color="black") +
  geom_hline(yintercept = nGene.low.cutoff.sd, color="black")

## Percent mito / nUMI plot
ggplot(pilot3.meta, aes(nUMI, percent.mito)) + 
  geom_hex(bins=100) + 
  scale_fill_distiller(palette = "Spectral", name="Cell frequencies") + 
  ylab("Fraction mitochondrial genes") + xlab("Number of reads") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=18)) + 
  geom_hline(yintercept = 0.05, colour="red")
mad