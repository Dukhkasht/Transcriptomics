rm(list =ls())
library(GEOquery)
library(affy)
library(tidyverse)
getwd()
#saved in gene analysis upgrade folder
untar("gse20966_RAW.tar",exdir= 'data3/')



# 1.A creating expression matrix using Raw file by affy
raw.data <- ReadAffy(celfile.path = "data3/")
normalised.data <- rma(raw.data)
gse20966 <- as.data.frame(exprs(normalised.data))



# 1.B creating expression matrix using series matrix

eSet <- getGEO("gse20966",
               destdir = '.',
               getGPL = F) 
gse20966 <- as.data.frame(exprs(eSet[[1]]))


#phenodata
gse <- getGEO("gse20966", GSEMatrix = TRUE)
pheno <- pData(phenoData(gse[[1]]))


# 2.A feature data to add gene symbol
feature.data <- gse$GSE20966_series_matrix.txt.gz@featureData@data
colnames(feature.data)
feature.data <- feature.data[,c(1,12)]


gse20966 <- gse20966 %>%
  rownames_to_column(var = 'ID')%>%
  inner_join(., feature.data,by = 'ID')
colnames(gse20966)
colnames(gse20966) <- gsub('\\..*','', colnames(gse20966))
gse20966=gse20966[!duplicated(gse20966$`Gene Symbol`),]
rownames(gse20966)= gse20966$`Gene Symbol`
gse20966= gse20966[,-c(1)]

save(gse20966,file = "diabetes2.Rdata")




# 2.B To add gen symbol from annotation file
gse20966$ID = rownames(gse20966)
data1<- read.csv('id35.csv',sep=',',header=T)
db3 = merge(gse20966,data1,by.x='ID',by.y='ID')
colnames(db3)
db3 = db3[!duplicated(db3$Gene.symbol),]
rownames(db3) = db3$Gene.symbol
db3 = db3[,-c(1)]
save(db3,file = 'Expression38642.Rdata')