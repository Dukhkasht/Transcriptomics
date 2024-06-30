rm(list =ls())


library(dplyr)
library(sva)
library(bladderbatch)

#loading Files having expression dataa
load('Expression45887.Rdata')
load('Expression6872.Rdata')


#Renaming column names
colnames(gse6872)
colnames(gse6872)[colnames(gse6872) == "symbol"] <- "Gene Symbol"


#merging all expression data
merge_eset<-inner_join(gse6872,gse45887, by = "Gene Symbol")
rownames(merge_eset)=merge_eset$'Gene Symbol'
print(colnames(merge_eset))
merge_eset<-merge_eset[,-c(22)]


#changing column names 
names(merge_eset) <- sub("\\..*", "", names(merge_eset))

colnames(merge_eset)


#creating file for batch correction 
exp=as.matrix(merge_eset)

data1<-read.table("mfbatch1.csv",sep=',',header=T)

data1$X<- gsub('\\_.*','',data1$X)

modcombat = model.matrix(~1, data = data1)

batch = data1$batch
combat_edata = ComBat(dat=merge_eset, batch=batch, mod=modcombat,par.prior=TRUE, prior.plots=FALSE)
combat <- as.data.frame(combat_edata)

write.table(combat, "mfallcounts.txt", sep = "\t", quote = F)


#visualisation
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
