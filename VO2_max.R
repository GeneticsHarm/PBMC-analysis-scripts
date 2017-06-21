library(Rtsne)

load(file = "~/Downloads/voor harm.RData")
all.real <- as.matrix(all[,c("Weight", "Length", "VO2max", "BMI", "RER", "HRmax", "WRmax")])
head(all.real)

all
all.real.pca$scores[5,]
all.real.pca$loadings

plot(data.a, data.b, col=factor(data.c))
plot(all.real.pca$scores[,1], all.real.pca$scores[,2], pch=20, cex=0.5, col = factor(all$Organisation))
plot(all.real.pca$scores[,1], all.real.pca$scores[,2], pch=20, cex=0.5, col = factor(all$Gender))
legend("topleft", legend=levels(factor(all$Organisation)), text.col=seq_along(levels(factor(all$Organisation))))

tsne_out <- Rtsne(all.real, perplexity=5, max_iter=500, check_duplicates=F, verbose=T)
qplot(tsne_out$Y[,1],tsne_out$Y[,2],data=as.data.frame(Ma))+labs(list(title="Barnes-Hut t-SNE",x="tSNE_1",y="tSNE_2",colour="condition"))+theme(legend.position="none")


all.real.pca$loadings
lines