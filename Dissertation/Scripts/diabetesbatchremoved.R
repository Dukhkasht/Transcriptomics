getwd()
load("expression20966.Rdata")
load("expression25724.Rdata")
library(dplyr)
library(sva)
library(bladderbatch)

GSE25724$`Gene Symbol` <- str_split(GSE25724$`Gene Symbol`, " /// ", simplify = TRUE)[, 1]
GSE25724=GSE25724[!duplicated(GSE25724$`Gene Symbol`),]
rownames(GSE25724)= GSE25724$`Gene Symbol`
#renaming
colnames(GSE25724)
colnames(gse38642)[colnames(gse38642) == "Gene.symbol"] <- "Gene Symbol"


#merging files on basis of gene symbol column
merge_eset<-inner_join(GSE25724,gse20966,by="Gene Symbol")

rownames(merge_eset)=merge_eset$'Gene Symbol'
print(colnames(merge_eset))
merge_eset<-merge_eset[,-c(14)]


exp=as.matrix(merge_eset)

data1<-read.table("diabetesnewbatch.csv",sep=',',header=T)
modcombat = model.matrix(~1, data = data1)
batch = data1$batch
combat_edata = ComBat(dat=merge_eset, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
write.table(combat_edata, "diabetesagaincounts.txt", sep = "\t", quote = F)

n <- as.data.frame(t(merge_eset))
name <- row.names(n)
par(las = 2, cex.lab = 0.5, cex.axis = 0.45, col = "black")

par(mar=c(6,2,1,1))
boxplot(as.data.frame(merge_eset),cex=0.2,las=2,col="red",main="Original",labels = name)
boxplot(as.data.frame(combat_edata),cex=0.2,las=2,col="red",main="Batch corrected",labels=name)
library("FactoMineR")
library("factoextra")
pca_plot = function(dddd,ggggg){
  library("FactoMineR")
  library("factoextra")
  df.pca <- PCA(t(dddd), graph = FALSE)
  fviz_pca_ind(df.pca,
               #axes = c(2,3),
               geom.ind = "point",
               col.ind = ggggg ,
               addEllipses = TRUE,
               legend.title="Groups"
  )
}

pca_plot(combat_edata,factor(data1$GSE))
pca_plot(merge_eset,factor(data1$GSE))

