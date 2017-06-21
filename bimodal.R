

expressionMeans <- rowMeans(as.matrix(combined.major@raw.data))
exp.filtered <- expressionMeans[expressionMeans>0.1]
length(exp.filtered)
exp.filtered

expressionMeans["MS4A1"]

par(mfrow=c(1,1))

hist(combined.major@data["MALAT1",], breaks = 100, main = "MALAT1 ~700 reads/cell (>100 = 22 genes)",  ylim = c(0,100))
hist(combined.major@data["RPL5",], breaks = 100, main = "RPL5 ~50 reads/cell (60 genes)")
hist(combined.major@data["LIMD2",], breaks = 100, main = "LIMD2 ~10 reads/cell (203 genes)", ylim = c(0,2000))
hist(combined.major@data["CSTB",], breaks = 100, main = "CSTB ~5 reads/cell (368 genes)", ylim = c(0,1000))
hist(combined.major@data["MS4A1",], breaks = 100, main = "MS4A1 ~2 reads/cell", ylim = c(0,300))
hist(combined.major@data["GNLY",], breaks = 100, main = "GNLY", ylim = c(0,500), xlab = "Log expression")
hist(combined.major@data["IL2RB",], breaks = 100, main = "IL2RB ~0.5 reads/cell (3194 genes)", ylim = c(0,100))
hist(combined.major@data["CD3E",], breaks = 100, main = "ADARB1 ~0.1 reads/cell (7796 genes)", ylim = c(0,900))
