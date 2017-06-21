eqtls.deepsage <- read.table("~/Documents/Pilot 3/eQTL/DeepSAGE_results.csv", sep = ";", header = T, dec = ",")
annotation <- read.table("~/Documents/Pilot 3/eQTL/singleCell-annotation-stripped.tsv", header=T, sep="\t", fill=T)

eqtls.deepsage[1:10,]
annotation.matched <- annotation[match(eqtls.deepsage$HGNCName, annotation$Symbol),]
annotation.matched$Symbol[1:10]
eqtls.deepsage$HGNCName[1:10]

snp.probe <- cbind.data.frame(eqtls.deepsage$SNPName, annotation.matched$Ensembl)
write.table(snp.probe, file="~/Desktop/DeepSAGE_probesnpfile_ensebl.tsv", sep = "\t", quote = F, row.names = F, col.names = F)
